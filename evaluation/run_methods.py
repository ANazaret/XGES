import inspect
import os
import shutil
import subprocess
import sys
import time
import uuid
from tempfile import gettempdir

import networkx as nx
import numpy as np
import pandas as pd
import sempler

from bic_score import BICScorer
from pdag import compute_ecshd_from_graph, PDAG

DEBUG = False


# like sdcd: benchmark suit
# > different simulations
# > different algorithms
# later:
# > intervention


def run_ges(data, alpha=1.0, phases=("forward", "backward", "turning")):
    from cdt.causality.graph import GES

    ges = GES()
    data = pd.DataFrame(data)
    graph = ges.predict(data, phases=phases, alpha=alpha)
    return graph


def run_fges(data, alpha):
    # save cwd
    cwd = os.getcwd()
    # required for java stuff
    os.chdir("py_tetrad/pytetrad")
    # required for python imports
    sys.path.append(os.path.join(os.path.dirname(__file__), "py_tetrad/pytetrad"))

    import jpype.imports

    try:
        jpype.startJVM(classpath=[f"resources/tetrad-current.jar"])
    except OSError:
        pass
    from py_tetrad.pytetrad.tools import TetradSearch as ts

    df = pd.DataFrame(data)
    search = ts.TetradSearch(df)
    # 1 = BIC; 2 = modified BIC from (Nandy, Hauser, Maathuis 2018)
    search.use_sem_bic(penalty_discount=alpha, sem_bic_rule=1)
    ## Run the search
    search.run_fges(parallelized=True)
    # restore cwd
    os.chdir(cwd)

    # we encode the tails of arrows as 1
    # undirected edges have two tails, so they are correctly encoded (1,1)
    # then we transpose the matrix
    graph_adjacency_matrix = search.get_graph_to_matrix(nullEpt=0, circleEpt=0, arrowEpt=0, tailEpt=1)

    return nx.DiGraph(graph_adjacency_matrix.values.T)


def run_xges(data, alpha, intervention_labels=None, threshold=True, ground_truth=None, extra_optim=0):
    # 1) save the data to a tmp location, generate a unique tmp file name (e.g. using uuid)

    if DEBUG:
        tmp_path = os.path.join("./", "xges_" + str(uuid.uuid4()))
    else:
        tmp_path = os.path.join(gettempdir(), "xges_" + str(uuid.uuid4()))

    os.makedirs(tmp_path, exist_ok=True)
    in_file = os.path.join(tmp_path, "data.npy")
    intervention_labels_file = os.path.join(tmp_path, "intervention_labels.npy")
    out_file = os.path.join(tmp_path, "result.csv")
    stats_file = os.path.join(tmp_path, "stats.csv")
    truth_file = os.path.join(tmp_path, "truth.csv")
    if isinstance(data, pd.DataFrame):
        data = data.values
    data = data.astype(np.float64)
    np.save(in_file, data)
    if intervention_labels is not None:
        intervention_labels = intervention_labels.astype(np.int32)
        np.save(intervention_labels_file, intervention_labels)

    # 2) run xges by calling the c++ binary: xges -f <tmp location> -a <alpha>
    xges_path = "../cmake-build-release/ges_cpp"
    cmd = [xges_path, "--input", in_file, "-a", str(alpha), "--output", out_file, "--stats", stats_file]
    if intervention_labels is not None:
        cmd += ["--interventions", intervention_labels_file]
    if threshold:
        cmd += ["-t"]
    if ground_truth is not None:
        df = pd.DataFrame(nx.to_numpy_array(ground_truth, nodelist=list(range(data.shape[1]))))
        df.to_csv(truth_file, index=False)
        cmd += ["-g", truth_file]

    if extra_optim:
        cmd += ["-o", str(extra_optim)]

    if DEBUG:
        print(" ".join(cmd))
        exit()
    start_time = time.time()
    print(" ".join(cmd))
    error_code = subprocess.run(cmd, timeout=60 * 60 * 24)
    if error_code.returncode != 0:
        print("error running xges")
        # show the error message, the cmd used and the tmp location
        print(error_code)
        print(" ".join(cmd))
        print(tmp_path)
        raise RuntimeError("error running xges")
    xges_time = time.time() - start_time

    # 3) load the result from the tmp location, pd.read_csv(out_file)
    graph = nx.DiGraph(pd.read_csv(out_file, header=0).values)
    xges_stats = pd.read_csv(stats_file, index_col=0, header=None).T.to_dict("records")[0]
    print(xges_stats)
    xges_stats["time_python"] = xges_time

    # 4) remove the tmp location
    if not DEBUG:
        shutil.rmtree(tmp_path)

    return graph, xges_stats


def compare_methods_given_data(data, true_graph, alpha=1.0, optim_xges=(0, 1, 2), intervention_labels=None):
    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data)
    # run ges, xges, pc, fges, ...

    res = dict()

    # # 1) run ges cdt (original R)
    # tmp_start_time = time.time()
    # ges_r_graph = run_ges(data)
    # ges_r_time = time.time() - tmp_start_time
    # ges_r_ecshd = compute_ecshd_from_graph(true_graph, ges_r_graph)
    # ges_r_shd = PDAG.from_digraph(ges_r_graph).shd_against_dag(true_graph)
    # res["ecshd_ges_r"] = ges_r_ecshd
    # res["time_ges_r"] = ges_r_time
    # res["shd_ges_r"] = ges_r_shd

    best_score = -np.inf
    associated_shd = None
    associated_ecshd = None
    # 2a) run xges
    tmp_start_time = time.time()
    xges_graph, stats = run_xges(data, alpha, intervention_labels=intervention_labels)
    xges_time = time.time() - tmp_start_time
    xges_ecshd = compute_ecshd_from_graph(true_graph, xges_graph)
    xges_shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
    res["ecshd_xges"] = xges_ecshd
    res["time_xges"] = xges_time
    res["shd_xges"] = xges_shd
    res["nedges_xges"] = len(xges_graph.edges)
    if stats["score"] > best_score:
        best_score = stats["score"]
        associated_shd = xges_shd
        associated_ecshd = xges_ecshd

    # # 2b) run xges with threshold
    # tmp_start_time = time.time()
    # xges_graph, stats = run_xges(data, alpha, threshold=True)
    # xges_time = time.time() - tmp_start_time
    # xges_ecshd = compute_ecshd_from_graph(true_graph, xges_graph)
    # xges_shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
    # res["ecshd_xges_t"] = xges_ecshd
    # res["time_xges_t"] = xges_time
    # res["shd_xges_t"] = xges_shd
    # if stats["score"] > best_score:
    #     best_score = stats["score"]
    #     associated_shd = xges_shd
    #     associated_ecshd = xges_ecshd
    #
    # # 2c) run xges with threshold and extra optim 1
    # tmp_start_time = time.time()
    # xges_graph, stats = run_xges(data, alpha, threshold=True, extra_optim=1)
    # xges_time = time.time() - tmp_start_time
    # xges_ecshd = compute_ecshd_from_graph(true_graph, xges_graph)
    # xges_shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
    # res["ecshd_xges_t_o1"] = xges_ecshd
    # res["time_xges_t_o1"] = xges_time
    # res["shd_xges_t_o1"] = xges_shd
    # if stats["score"] > best_score:
    #     best_score = stats["score"]
    #     associated_shd = xges_shd
    #     associated_ecshd = xges_ecshd

    for o in optim_xges:
        # 2d) run xges with threshold and extra optim 2
        tmp_start_time = time.time()
        xges_graph, stats = run_xges(
            data, alpha, threshold=True, extra_optim=o, intervention_labels=intervention_labels
        )
        xges_time = time.time() - tmp_start_time
        xges_ecshd = compute_ecshd_from_graph(true_graph, xges_graph)
        xges_shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
        print("SHD", xges_shd, "EC-SHD", xges_ecshd, "TIME", xges_time, "SCORE", stats["score"])
        res[f"ecshd_xges_t_o{o}"] = xges_ecshd
        res[f"time_xges_t_o{o}"] = xges_time
        res[f"shd_xges_t_o{o}"] = xges_shd
        res[f"nedges_xges_t_o{o}"] = len(xges_graph.edges)
        if stats["score"] > best_score:
            best_score = stats["score"]
            associated_shd = xges_shd
            associated_ecshd = xges_ecshd

    # # 2e) run xges no threshold and extra optim 1
    # tmp_start_time = time.time()
    # xges_graph, stats = run_xges(data, alpha, threshold=False, extra_optim=1)
    # xges_time = time.time() - tmp_start_time
    # xges_ecshd = compute_ecshd_from_graph(true_graph, xges_graph)
    # xges_shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
    # res["ecshd_xges_o1"] = xges_ecshd
    # res["time_xges_o1"] = xges_time
    # res["shd_xges_o1"] = xges_shd
    # if stats["score"] > best_score:
    #     best_score = stats["score"]
    #     associated_shd = xges_shd
    #     associated_ecshd = xges_ecshd

    # 2z) pick the best score
    res["shd_xges_best"] = associated_shd
    res["ecshd_xges_best"] = associated_ecshd

    # 3) run fges
    # tmp_start_time = time.time()
    # fges_graph = run_fges(data, alpha)
    # fges_time = time.time() - tmp_start_time
    # fges_ecshd = compute_ecshd_from_graph(true_graph, fges_graph)
    # fges_shd = PDAG.from_digraph(fges_graph).shd_against_dag(true_graph)
    # res["ecshd_fges"] = fges_ecshd
    # res["time_fges"] = fges_time
    # res["shd_fges"] = fges_shd

    return res


def run_method_given_data(method_name, data, true_graph, **kwargs):
    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data)

    method_func = {
        "ges": run_ges,
        "ges+turn": run_ges,
        "fges": run_fges,
        "xgesZ": run_xges,
        "xges1Z": run_xges,
        "xges0": run_xges,
        "xges1": run_xges,
        "xges2": run_xges,
    }[method_name]

    if method_name == "ges":
        kwargs["phases"] = ("forward", "backward")
    elif method_name == "ges+turn":
        kwargs["phases"] = ("forward", "backward", "turning")

    if method_name == "xgesZ":
        kwargs["threshold"] = False
        kwargs["extra_optim"] = 0
    if method_name == "xges1Z":
        kwargs["threshold"] = False
        kwargs["extra_optim"] = 1
    if method_name == "xges0":
        kwargs["threshold"] = True
        kwargs["extra_optim"] = 0
    elif method_name == "xges1":
        kwargs["threshold"] = True
        kwargs["extra_optim"] = 1
    elif method_name == "xges2":
        kwargs["threshold"] = True
        kwargs["extra_optim"] = 2

    tmp_start_time = time.time()
    graph = run_filter_kwargs(method_func, kwargs, data=data, true_graph=true_graph)
    res = dict()
    if isinstance(graph, tuple):
        graph, stats = graph
        res["score"] = stats["score"]
    res["time"] = time.time() - tmp_start_time

    # compute metrics
    res["ecshd"] = compute_ecshd_from_graph(true_graph, graph)
    res["shd"] = PDAG.from_digraph(graph).shd_against_dag(true_graph)
    res["n_edges_predicted"] = len(graph.edges)
    res["n_edges_true"] = len(true_graph.edges)
    # score analysis
    bs = BICScorer(
        data.shape[0],
        data.shape[1],
        alpha=kwargs["alpha"],
        data=data,
        # intervention_labels=intervention_labels,
    )
    res["score_true_graph"] = bs.score_of_dag(true_graph)
    res["score_check"] = bs.score_of_dag(graph)
    precision = 0  # in the predicted graph, proportion of correct edges
    recall = 0  # among all true edges, proportion of those who are predicted
    n_predicted_edges = 0
    for edge in graph.edges:
        if graph.has_edge(edge[1], edge[0]):
            if edge[0] < edge[1]:
                n_predicted_edges += 1
        else:
            n_predicted_edges += 1
        if true_graph.has_edge(*edge):
            precision += 1
    precision /= n_predicted_edges

    n_true_edges = len(true_graph.edges)
    for edge in true_graph.edges:
        if graph.has_edge(*edge):
            recall += 1
    recall /= n_true_edges

    res["precision"] = precision
    res["recall"] = recall
    res["f1"] = 2 * precision * recall / (precision + recall)
    return res


def compare_methods(
    n_seeds=5,
    d_values=(10, 20, 50, 100),
    edges_per_d_values=(1, 2, 3, 4),
    alphas=(1.0, 2.0),
    optim_xges=(0, 1, 2),
    n_inter=(0, 100),
):
    scores = []
    for d in d_values:
        print("d", d)
        for edges_per_d in edges_per_d_values:
            print("d", d, "edges_per_d", edges_per_d)
            for seed in range(n_seeds):
                for n_inter_ in n_inter:
                    print("d", d, "edges_per_d", edges_per_d, "seed", seed)
                    data, true_graph, intervention_labels = simulation1(
                        d,
                        n_control=10_000,
                        n_inter=n_inter_,
                        edges_per_d=edges_per_d,
                        seed=seed,
                        normalize=True,
                    )
                    # if n_inter_ > 0:
                    #     # one hot encode the interventions
                    #     tmp = np.eye(d + 1)[intervention_labels]
                    #     tmp = tmp[:, :-1]
                    #     print(np.cov(np.concatenate([data, tmp], axis=1), rowvar=False, bias=True))

                    for alpha in alphas:
                        print("d", d, "edges_per_d", edges_per_d, "seed", seed, "alpha", alpha)
                        res = compare_methods_given_data(
                            data,
                            true_graph,
                            alpha=alpha,
                            optim_xges=optim_xges,
                            intervention_labels=intervention_labels,
                        )
                        res["d"] = d
                        res["seed"] = seed
                        res["edges_per_d"] = edges_per_d
                        res["alpha"] = alpha
                        res["n_inter"] = n_inter_
                        res["n_obs"] = data.shape[0]
                        scores.append(res)

    scores = pd.DataFrame(scores)
    # sort columns
    scores = scores[list(sorted(scores.columns))]
    print(scores.groupby(["d", "edges_per_d"]).mean().to_string())
    print(scores.groupby(["d", "edges_per_d", "alpha"]).mean().to_string())
    print(scores.groupby(["d", "edges_per_d", "alpha", "n_inter"]).mean().to_string())

    # print(scores.groupby(["d", "edges_per_d"]).median().to_string())
    return scores


def run_filter_kwargs(func, all_kwargs, **kwargs):
    func_args = inspect.signature(func).parameters
    all_kwargs = all_kwargs.copy()
    all_kwargs.update(kwargs)
    return func(**{k: v for k, v in all_kwargs.items() if k in func_args})


def worker(kwargs_):
    data, true_graph, intervention_labels = run_filter_kwargs(simulation1, kwargs_)
    res = kwargs_.copy()
    kwargs_["intervention_labels"] = intervention_labels
    res_2 = run_method_given_data(data=data, true_graph=true_graph, **kwargs_)
    res.update(res_2)
    res["n_obs"] = data.shape[0]
    print(res["shd"])
    print(res)
    return res


def compare_methodsv2(n_jobs=1, **kwargs):
    scores = []
    # do all combinations of the parameters
    from itertools import product
    import multiprocessing

    if n_jobs == 1:
        scores = [worker(dict(zip(kwargs.keys(), kwargs_))) for kwargs_ in product(*kwargs.values())]
    else:
        with multiprocessing.Pool(n_jobs) as p:
            scores = p.map(
                worker, [dict(zip(kwargs.keys(), kwargs_)) for kwargs_ in product(*kwargs.values())]
            )

    # for kwargs_ in product(*kwargs.values()):
    #     kwargs_ = dict(zip(kwargs.keys(), kwargs_))
    #     print(kwargs_)
    #     data, true_graph, intervention_labels = run_filter_kwargs(simulation1, kwargs_)
    #     res = kwargs_.copy()
    #     kwargs_["intervention_labels"] = intervention_labels
    #     res_2 = run_method_given_data(data=data, true_graph=true_graph, **kwargs_)
    #     res.update(res_2)
    #     res["n_obs"] = data.shape[0]
    #     print(res["shd"])
    #     scores.append(res)
    scores = pd.DataFrame(scores)
    # sort columns
    scores = scores[list(sorted(scores.columns))]
    return scores


def benchmark_xges_alpha_given_data(data, true_graph):
    res_stats = []
    for alpha in [0.1, 0.5, 1, 2, 5, 10]:
        graph, stats = run_xges(data, alpha)
        res_stats.append(
            {
                "alpha": alpha,
                "time": stats["time"],
                "ec-shd": compute_ecshd_from_graph(graph, true_graph),
                "shd": PDAG.from_digraph(graph).shd_against_dag(true_graph),
                "score": stats["score"],
                "score_increase": stats["score_increase"],
            }
        )

    return pd.DataFrame(res_stats)


def benchmark_xges_alpha():
    scores = []
    for d in [10, 20, 50, 100]:
        for seed in range(2):
            data, true_graph = simulation1(
                d,
                n_control=10000,
                n_inter=0,
                edges_per_d=int(d**0.5),
                seed=seed,
                normalize=True,
                intervention=False,
            )
            res = benchmark_xges_alpha_given_data(data, true_graph)
            res["d"] = d
            res["seed"] = seed
            scores.append(res)

    scores = pd.concat(scores)
    scores.to_csv("scores-xges-alphas.csv", index=False)


def experiment_more_observational():
    scores = []
    n_inter = 100
    n_obs_vals = [100, 1000, 10000, 100000]
    for d in [30]:
        for seed in range(10):
            xges_graphs = dict()
            df = []

            for n_obs in n_obs_vals:
                data, true_graph, intervention_labels = simulation1(
                    d,
                    n_control=n_obs,
                    n_inter=n_inter,
                    edges_per_d=4,
                    seed=seed,
                    normalize=True,
                    intervention=True,
                )
                alpha = int(np.sqrt(data.shape[0]))

                xges_graph, xges_stats = run_xges(data, alpha, intervention_labels)
                xges_graphs[n_obs] = xges_graph
            xges_graphs["true"] = true_graph

            for graph_idx in xges_graphs:
                # compute the score of each graph under each n_obs
                df_record = dict()
                df_record["graph"] = graph_idx

                for n_obs in n_obs_vals:
                    data, true_graph, intervention_labels = simulation1(
                        d,
                        n_control=n_obs,
                        n_inter=n_inter,
                        edges_per_d=4,
                        seed=seed,
                        normalize=True,
                        intervention=True,
                    )
                    bs = BICScorer(data.shape[0], data.shape[1], alpha=alpha, data=data, ignored_variables=[])
                    df_record[f"score_n_obs{n_obs}"] = bs.score_of_dag(xges_graphs[graph_idx])

                df_record["shd"] = PDAG.from_digraph(xges_graphs[graph_idx]).shd_against_dag(true_graph)
                df_record["ec-shd"] = compute_ecshd_from_graph(xges_graphs[graph_idx], true_graph)
                df_record["n_edges"] = len(xges_graphs[graph_idx].edges)
                df_record["alpha"] = alpha
                df_record["d"] = d
                df_record["seed"] = seed

                df.append(df_record)
            df = pd.DataFrame(df)
            scores.append(df)

    scores = pd.concat(scores)
    scores.to_csv("scores-xges-range-obs.csv", index=False)


def test_score_match():
    data, true_graph = simulation1(20, n_control=10000, n_inter=0, edges_per_d=3, seed=0, normalize=True)
    alpha = 1.0
    xges_graph, xges_stats = run_xges(data, alpha, threshold=False)

    bs = BICScorer(data.shape[0], data.shape[1], alpha=alpha, data=data, ignored_variables=[])
    cpp_score = xges_stats["score"]
    python_score = bs.score_of_dag(xges_graph)

    print(
        "python score",
        python_score,
        "cpp empty score",
        cpp_score,
        "shd",
        PDAG.from_digraph(xges_graph).shd_against_dag(true_graph),
    )
    print(true_graph.edges)
    print(list(nx.topological_sort(true_graph)))
    # plot graph in topological order

    # nx.draw_networkx(true_graph, with_labels=True)

    # plt.show()
    # data, true_graph, interventions_index = simulation1(
    #     20, n_control=10000, n_inter=100, edges_per_d=2, seed=0, normalize=True, intervention=True
    # )
    # alpha = 2
    # xges_graph, xges_stats = run_xges(data, alpha, interventions_index)
    #
    # bs = BICScorer(data.shape[0], data.shape[1], alpha=alpha, data=data, ignored_variables=[])
    # our_score = bs.score_of_dag(xges_graph)
    # true_graph_our_score = bs.score_of_dag(true_graph)
    # print(
    #     "our score",
    #     our_score,
    #     "xges score",
    #     xges_stats["score"],
    #     "true graph score",
    #     true_graph_our_score,
    # )


def test_threshold():
    res = []
    for d in [20, 30, 40, 50, 100]:
        for edges_per_d in [1, 2, 3, 4]:
            for seed in range(5):
                for n_obs in [1000, 10000, 100000]:
                    data, true_graph = simulation1(
                        # data, true_graph, ii = simulation1(
                        d,
                        n_control=n_obs,
                        n_inter=100,
                        edges_per_d=edges_per_d,
                        seed=seed,
                        normalize=True,
                        # intervention=True,
                    )
                    alpha = 1.0

                    best_score = -np.inf
                    associated_shd = None
                    for t in [True, False]:
                        # xges_graph, xges_stats = run_xges(data, alpha, ii, threshold=t)
                        xges_graph, xges_stats = run_xges(data, alpha, threshold=t)
                        cpp_score = xges_stats["score"]
                        shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
                        res.append(
                            {
                                "d": d,
                                "edges_per_d": edges_per_d,
                                "seed": seed,
                                "n_obs": n_obs,
                                "cpp_score": cpp_score,
                                "shd": shd,
                                "threshold": t,
                            }
                        )
                        if cpp_score > best_score:
                            best_score = cpp_score
                            associated_shd = shd
                    res.append(
                        {
                            "d": d,
                            "edges_per_d": edges_per_d,
                            "seed": seed,
                            "n_obs": n_obs,
                            "cpp_score": best_score,
                            "shd": associated_shd,
                            "threshold": "best",
                        }
                    )
    res = pd.DataFrame(res)
    print(res.groupby(["d", "n_obs", "threshold"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["d", "threshold"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["n_obs", "threshold"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["edges_per_d", "threshold"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["threshold"])[["cpp_score", "shd"]].mean().to_string())

    res.to_csv("scores-xges-threshold.csv", index=False)
    # res.to_csv("scores-xges-threshold-interventions.csv", index=False)

    # plot graph in topological order


def test_optim():
    res = []
    for d in [20, 30, 40, 50, 100]:
        for edges_per_d in [1, 2, 3, 4]:
            for seed in range(3):
                for n_obs in [1000, 10000, 100000]:
                    data, true_graph = simulation1(
                        # data, true_graph, ii = simulation1(
                        d,
                        n_control=n_obs,
                        n_inter=100,
                        edges_per_d=edges_per_d,
                        seed=seed,
                        normalize=True,
                        # intervention=True,
                    )
                    alpha = 1.0

                    for optim in [True, False]:
                        # xges_graph, xges_stats = run_xges(data, alpha, ii, threshold=t)
                        xges_graph, xges_stats = run_xges(data, alpha, threshold=True, extra_optim=optim)
                        cpp_score = xges_stats["score"]
                        shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
                        res.append(
                            {
                                "d": d,
                                "edges_per_d": edges_per_d,
                                "seed": seed,
                                "n_obs": n_obs,
                                "cpp_score": cpp_score,
                                "shd": shd,
                                "extra_optim": optim,
                            }
                        )
    res = pd.DataFrame(res)
    # print(res.groupby(["d", "n_obs", "extra_optim"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["d", "extra_optim"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["n_obs", "extra_optim"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["edges_per_d", "extra_optim"])[["cpp_score", "shd"]].mean().to_string())
    print(res.groupby(["extra_optim"])[["cpp_score", "shd"]].mean().to_string())

    res.to_csv("scores-xges-extra_optim.csv", index=False)


def test_reproducibility():
    data, true_graph = simulation1(
        50,
        n_control=10_000,
        n_inter=100,
        edges_per_d=2,
        seed=0,
        normalize=True,
        # intervention=True,
    )
    alpha = 1.0

    xges_graph, xges_stats = run_xges(data, alpha, ground_truth=true_graph, threshold=True, extra_optim=True)
    cpp_score = xges_stats["score"]
    shd = PDAG.from_digraph(xges_graph).shd_against_dag(true_graph)
    print(cpp_score, shd)


def simulation1(
    d,
    n_control=1000,
    n_inter=0,
    edges_per_d=2,
    seed=0,
    normalize=True,
    max_variance=1,
    shift_value=-1,
    positive=True,
):
    rng = np.random.RandomState(seed)
    graph = nx.gnp_random_graph(d, edges_per_d / (d - 1), directed=True, seed=seed)
    # make it acyclic
    true_dag = nx.DiGraph()
    true_dag.add_nodes_from(graph.nodes())
    true_dag.add_edges_from([(u, v) if u < v else (v, u) for (u, v) in graph.edges() if u != v])
    # shuffle the nodes
    permutation = rng.permutation(d)
    true_dag = nx.relabel_nodes(true_dag, {i: permutation[i] for i in range(d)})

    true_adjacency = nx.to_numpy_array(true_dag, nodelist=list(range(d)))
    weights = true_adjacency * rng.uniform(
        1, 3, true_adjacency.shape
    )  # * rng.choice([-1, 1], true_adjacency.shape)
    if not positive:
        weights *= rng.choice([-1, 1], true_adjacency.shape)
    if normalize:
        # not sure if it is axis=0 or axis=1 ...
        # BE CAREFUL WITH the choice(-1, 1) above; sum would need to be abs or something
        weights = weights / (np.abs(weights).sum(axis=0, keepdims=True) + 1e-5)

    for e in true_dag.edges:
        true_dag.edges[e]["w"] = weights[e[0], e[1]]
    scm = sempler.LGANM(weights, (-2, 2), (0, max_variance), random_state=seed)
    data = scm.sample(
        n=n_control,
        random_state=seed,
    )

    data = data.astype(np.float64)
    if n_inter > 0:
        data = [data]
        intervention_labels = [-1] * n_control
        for i in range(d):
            intervention_labels += [i] * n_inter
            data.append(scm.sample(n=n_inter, random_state=i, shift_interventions={i: (shift_value, 0)}))
        data = np.concatenate(data, axis=0)
        intervention_labels = np.array(intervention_labels, dtype=np.int32)
        data = data.astype(np.float64)
        return data, true_dag, intervention_labels
    return data, true_dag, None


# DEBUG = True
if __name__ == "__main__":
    # data, _, _ = simulation1(20, n_control=10000, n_inter=0, edges_per_d=3, seed=0, normalize=True)
    # print(data.flags)
    # experiment_more_observational()
    # test_score_match()
    # benchmark_xges_alpha()
    # test_threshold()
    # test_optim()
    # test_reproducibility()

    # res = compare_methods(10, [30, 40, 50], [2, 3, 4], alphas=[2])
    # res.to_csv("res.csv", index=False)
    # exit()

    # res_inter = compare_methodsv2(
    #     seed=list(range(10)),
    #     n_control=[10_000],
    #     method_name=["xges0", "xges1"],
    #     d=[30, 40],
    #     edges_per_d=[3],
    #     alpha=[2],
    #     n_inter=[0, 10, 20, 30, 40, 50, 60, 70, 100, 150],
    # )
    # res_inter.to_csv("res-inter2.csv", index=False)

    # res_general = compare_methodsv2(
    #     seed=list(range(10)),
    #     n_control=[10_000],
    #     method_name=[
    #         "fges",
    #         "ges",
    #         "xges0",
    #         "xges1",
    #         "xges2",
    #     ],
    #     d=[20, 30, 40, 50, 60, 70],
    #     edges_per_d=[2, 3, 4],
    #     alpha=[2],
    #     n_inter=[0],
    # )
    # res_general.to_csv("res-general2.csv", index=False)
    # res_xges = compare_methodsv2(
    #     seed=list(range(10)),
    #     n_control=[10_000],
    #     method_name=[
    #         "xges0",
    #         "xges1",
    #     ],
    #     d=[50, 100, 200, 300, 400, 500, 1000, 2000],
    #     edges_per_d=[2, 3, 4],
    #     alpha=[2],
    #     n_inter=[0],
    # )
    # res_xges.to_csv("res-xges2.csv", index=False)

    # res_paper_0 = compare_methodsv2(
    #     n_jobs=8,
    #     seed=list(range(30)),
    #     n_control=[10_000],
    #     method_name=[
    #         "fges",
    #         "ges",
    #         "ges+turn",
    #         "xges0",
    #         "xges1",
    #         "xges2",
    #     ],
    #     d=[10, 25, 50, 75, 100, 125, 150, 175, 200],
    #     edges_per_d=[2, 3, 4],
    #     alpha=[2],
    #     n_inter=[0],
    #     max_variance=[0.5],
    # )
    # res_paper_0.to_csv("res-paper-1h.csv", index=False)

    # warm up the algorithms
    # compare_methodsv2(
    #     n_jobs=1,
    #     seed=list(range(1)),
    #     n_control=[10_000],
    #     method_name=[
    #         "fges",
    #         "ges",
    #         "xges0",
    #     ],
    #     # d=[10, 25, 50, 100],
    #     d=[10],
    #     edges_per_d=[1],
    #     alpha=[2],
    #     n_inter=[0],
    #     max_variance=[0.5],
    # )

    # res_paper_d = compare_methodsv2(
    #     n_jobs=1,
    #     seed=list(range(5)),
    #     n_control=[10_000],
    #     method_name=[
    #         # "fges",
    #         # "ges",
    #         # "ges+turn",
    #         "xges1Z",
    #         # "xges0",
    #         # "xges1",
    #         # "xges2",
    #     ],
    #     d=[10, 25, 50, 100, 200, 300],
    #     # d=[750],
    #     edges_per_d=[1, 2, 3, 4],
    #     alpha=[2],
    #     n_inter=[0],
    #     max_variance=[0.5],
    # )
    # res_paper_d.to_csv("res-paper-benchmark-p13.csv", index=False)
    # p1: d=[10, 25, 50, 100], edges_per_d=[1, 2, 3, 4] -> no warm up
    # p2: d=[200, 300], edges_per_d=[1, 2, 3, 4] -> no warm up
    # p3: d=[10, 25, 50, 100,], edges_per_d=[5] -> warm up
    # p4: d=[200, 300], edges_per_d=[5] -> warm up, no fges
    # p5: d=[500], edges_per_d=[1, 2, 3, 4, 5] -> warm up, no fges
    # p6: d=[750], edges_per_d=[1, 2, 3, 4, 5] -> warm up, no fges
    # p7: d=[500], edges_per_d=[1, 3] -> warm up, fges only
    # p8: d=[750], edges_per_d=[1, 3] -> warm up, fges only
    # p9: d=[10, 25, 50, 100,200,300,500,750], edges_per_d=[1, 3] -> warm up, fges only parallelized
    # p10: d=[10, 25, 50, 100,200,300,500,750], edges_per_d=[5] -> warm up, fges only parallelized
    # p11: d=[500, 750], edges_per_d=[2] -> warm up, fges only parallelized
    # p12: d=[10, 25, 50, 100, 200, 300], edges_per_d=[2,3,4] -> warm up xgesZ
    # p13: d=[10, 25, 50, 100, 200, 300], edges_per_d=[2,3,4] -> warm up xges1Z

    # res_paper_1 = compare_methodsv2(
    #     n_jobs=8,
    #     seed=list(range(5)),
    #     n_control=[10_000],
    #     method_name=[
    #         "fges",
    #         "ges",
    #         "ges+turn",
    #         "xges0",
    #         "xges1",
    #         "xges2",
    #     ],
    #     d=[10, 50, 100],
    #     edges_per_d=[2, 3, 4],
    #     alpha=[2],
    #     n_inter=[0],
    #     max_variance=[0.1, 0.2, 0.5, 1],
    #     positive=[True, False],
    # )
    # # 1e: vary alpha
    # # 1f: vary max_variance / positive
    # res_paper_1.to_csv("res-paper-1f.csv", index=False)

    # res_paper_1 = compare_methodsv2(
    #     n_jobs=8,
    #     seed=list(range(5)),
    #     n_control=[10_000],
    #     method_name=[
    #         "fges",
    #         "ges",
    #         "ges+turn",
    #         "xges0",
    #         "xges1",
    #         # "xges2",
    #     ],
    #     d=[25, 50],
    #     edges_per_d=[2, 3, 4],
    #     alpha=[2],
    #     n_inter=[0],
    #     max_variance=[0.5],
    #     positive=[True, False],
    # )
    # # 1f: vary positive
    # res_paper_1.to_csv("res-paper-positive.csv", index=False)

    # res_paper_2 = compare_methodsv2(
    #     n_jobs=8,
    #     seed=list(range(10)),
    #     n_control=[10_000],
    #     method_name=[
    #         "fges",
    #         "ges",
    #         "ges+turn",
    #         "xges0",
    #         "xges1",
    #         "xges2",
    #     ],
    #     d=[10, 50, 100],
    #     edges_per_d=[2, 3, 4],
    #     alpha=[1, 2, 5, 7, 10, 15, 20, 30],
    #     n_inter=[0],
    #     max_variance=[0.2],
    # )
    # res_paper_2.to_csv("res-paper-1e.csv", index=False)

    # res_paper_3 = compare_methodsv2(
    #     n_jobs=8,
    #     seed=list(range(10)),
    #     # n_control=[100, 500, 1_000, 5_000, 10_000, 50_000, 100_000],
    #     n_control=[500_000, 1_000_000, 5_000_000],
    #     method_name=[
    #         # "fges",
    #         # "ges",
    #         # "ges+turn",
    #         "xges0",
    #         "xges1",
    #         "xges2",
    #     ],
    #     d=[50],
    #     edges_per_d=[3],
    #     alpha=[2],
    #     n_inter=[0],
    #     max_variance=[0.5],
    # )
    # res_paper_3.to_csv("res-paper-1gb.csv", index=False)

    tmp = pd.read_csv("../synthetic_data.csv")
    g = run_xges(tmp, 2, threshold=False, extra_optim=1)
    print(g)
