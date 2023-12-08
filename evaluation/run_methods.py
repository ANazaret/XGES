import os
import shutil
import subprocess
import time
import uuid
from tempfile import gettempdir

import networkx as nx
import numpy as np
import pandas as pd
import sempler
from cdt.causality.graph import GES

from bic_score import BICScorer
from pdag import compute_ecshd_from_graph, PDAG


def simulation1(d, n_control=1000, n_inter=100, edges_per_d=10, seed=0, normalize=True, intervention=False):
    rng = np.random.RandomState(seed)
    graph = nx.gnp_random_graph(d, edges_per_d / d, directed=True, seed=seed)
    # make it acyclic
    true_dag = nx.DiGraph()
    true_dag.add_nodes_from(graph.nodes())
    true_dag.add_edges_from([(u, v) if u < v else (v, u) for (u, v) in graph.edges() if u != v])
    # shuffle the nodes
    permutation = rng.permutation(d)
    true_dag = nx.relabel_nodes(true_dag, {i: permutation[i] for i in range(d)})

    true_adjacency = nx.to_numpy_array(true_dag, nodelist=list(range(d)))
    weights = (
        true_adjacency * rng.uniform(1, 3, true_adjacency.shape) * rng.choice([-1, 1], true_adjacency.shape)
    )
    if normalize:
        # not sure if it is axis=0 or axis=1 ...
        weights = weights / (weights.sum(axis=0, keepdims=True) + 1e-5)

    for e in true_dag.edges:
        true_dag.edges[e]["w"] = weights[e[0], e[1]]
    scm = sempler.LGANM(weights, (-2, 2), (0, 1), random_state=seed)
    data = scm.sample(
        n=n_control,
        random_state=seed,
    )

    data = data.astype(np.float64)
    if intervention:
        data = [data]
        intervention_labels = [-1] * n_control
        for i in range(d):
            intervention_labels += [i] * n_inter
            data.append(scm.sample(n=n_inter, random_state=i, shift_interventions={i: (-0.1, 0)}))
        data = np.concatenate(data, axis=0)
        intervention_labels = np.array(intervention_labels, dtype=np.int32)
        data = data.astype(np.float64)
        return data, true_dag, intervention_labels
    return data, true_dag


# like sdcd: benchmark suit
# > different simulations
# > different algorithms
# later:
# > intervention


def run_ges(data, alpha, **kwargs):
    ges = GES()
    output = ges.predict(data, alpha=alpha, **kwargs)
    return output.graph


def run_xges(data, alpha, intervention_labels=None, threshold=False):
    # 1) save the data to a tmp location, generate a unique tmp file name (e.g. using uuid)

    tmp_path = os.path.join(gettempdir(), "xges_" + str(uuid.uuid4()))
    os.makedirs(tmp_path, exist_ok=True)
    in_file = os.path.join(tmp_path, "data.npy")
    intervention_labels_file = os.path.join(tmp_path, "intervention_labels.npy")
    out_file = os.path.join(tmp_path, "result.csv")
    stats_file = os.path.join(tmp_path, "stats.csv")
    data = data.astype(np.float64)
    np.save(in_file, data)
    if intervention_labels is not None:
        intervention_labels = intervention_labels.astype(np.int32)
        np.save(intervention_labels_file, intervention_labels)

    # 2) run xges by calling the c++ binary: xges -f <tmp location> -a <alpha>
    xges_path = "../cmake-build-release/ges_cpp"
    cmd = [xges_path, "-i", in_file, "-a", str(alpha), "-o", out_file, "--stats", stats_file]
    if intervention_labels is not None:
        cmd += ["--interventions", intervention_labels_file]
    if threshold:
        cmd += ["-t"]
    start_time = time.time()
    error_code = subprocess.run(cmd)
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
    shutil.rmtree(tmp_path)

    return graph, xges_stats


def compare_methods_given_data(data, true_graph, alpha=1):
    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data)
    # run ges, xges, pc, fges, ...

    # 1) run ges cdt (original R)
    ges_r = GES()
    tmp_start_time = time.time()
    ges_r_graph = ges_r.predict(data)
    ges_r_time = time.time() - tmp_start_time
    ges_r_ecshd = compute_ecshd_from_graph(true_graph, ges_r_graph)

    # 2) run xges
    tmp_start_time = time.time()
    xges_graph, _ = run_xges(data, alpha)
    xges_time = time.time() - tmp_start_time
    xges_ecshd = compute_ecshd_from_graph(true_graph, xges_graph)

    print()
    print("ges_r_shd", ges_r_ecshd)
    print("xges_shd", xges_ecshd)

    return {
        "ges_r_ecshd": ges_r_ecshd,
        "xges_ecshd": xges_ecshd,
        "ges_r_time": ges_r_time,
        "xges_time": xges_time,
    }


def compare_methods():
    scores = []
    for d in [10, 20, 50, 100, 200]:
        for seed in range(10):
            data, true_graph = simulation1(
                d, n_control=10000, n_inter=None, edges_per_d=4, seed=seed, normalize=True, intervention=False
            )
            res = compare_methods_given_data(data, true_graph)
            res["d"] = d
            res["seed"] = seed
            scores.append(res)

    scores = pd.DataFrame(scores)
    scores.to_csv("scores.csv", index=False)
    print(scores)
    print(scores.groupby("d").mean().to_string())


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
                    data, true_graph, ii = simulation1(
                        d,
                        n_control=n_obs,
                        n_inter=100,
                        edges_per_d=edges_per_d,
                        seed=seed,
                        normalize=True,
                        intervention=True,
                    )
                    alpha = 1.0

                    best_score = -np.inf
                    associated_shd = None
                    for t in [True, False]:
                        xges_graph, xges_stats = run_xges(data, alpha, ii, threshold=t)
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

    # res.to_csv("scores-xges-threshold.csv", index=False)
    res.to_csv("scores-xges-threshold-interventions.csv", index=False)

    # plot graph in topological order


if __name__ == "__main__":
    # experiment_more_observational()
    # test_score_match()
    # benchmark_xges_alpha()
    test_threshold()
