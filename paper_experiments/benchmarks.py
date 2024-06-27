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

from .bic_score import BICScorer
from .pdag import PDAG, compute_shd_cpdag


def run_ges(data, alpha=2.0, phases=("forward", "backward", "turning")):
    """Run the GES algorithm on the given data.

    This is a wrapper around the GES algorithm from the R package pcalg.

    Parameters
    ----------
    data : pd.DataFrame
        The data to run the algorithm on.
    alpha : float (default=2.0)
        The alpha parameter for the BIC score.
    phases : tuple (default=("forward", "backward", "turning"))
        The phases to run the algorithm with, corresponding respectively to: insertions, deletions, and reversals.

    Returns
    -------
    nx.DiGraph
        The graph learned by the algorithm. Undirected edges are represented as two directed edges.
    """
    from GES_custom import GES

    ges = GES()
    data = pd.DataFrame(data)
    graph = ges.predict(data, phases=phases, alpha=alpha)
    return graph


def run_fges(data, alpha=2.0):
    """Run the FGES algorithm on the given data.

    This is a wrapper around the FGES algorithm from the Java software Tetrad.

    Parameters
    ----------
    data : pd.DataFrame
        The data to run the algorithm on.
    alpha : float (default=2.0)
        The alpha parameter for the BIC score.

    Returns
    -------
    nx.DiGraph
        The graph learned by the algorithm. Undirected edges are represented as two directed edges.
    """
    # save cwd
    cwd = os.getcwd()
    # required for java stuff
    os.chdir("pytetrad")
    # required for python imports
    sys.path.append(os.path.join(os.path.dirname(__file__), "pytetrad"))

    import jpype.imports

    try:
        jpype.startJVM(classpath=["resources/tetrad-current.jar"])
    except OSError:
        pass
    from tools import TetradSearch as ts

    df = pd.DataFrame(data)
    search = ts.TetradSearch(df)
    # 1 = BIC; 2 = modified BIC from (Nandy, Hauser, Maathuis 2018)
    search.use_sem_bic(penalty_discount=alpha, sem_bic_rule=1)
    ## Run the search
    search.run_fges(parallelized=True)
    # restore cwd
    os.chdir(cwd)

    # encode the tails of arrows as 1
    # undirected edges have two tails, so they are encoded (1,1), corresponding to two directed edges
    # then we transpose the matrix to have the 1 toward the head of the arrow
    graph_adjacency_matrix = search.get_graph_to_matrix(
        nullEpt=0, circleEpt=0, arrowEpt=0, tailEpt=1
    )
    return nx.DiGraph(graph_adjacency_matrix.values.T)


def run_xges(data, alpha=2.0, extended_search=True, ground_truth=None, baseline="", verbose=0):
    """Run the XGES algorithm on the given data.

    This is a wrapper around the XGES algorithm from the C++ software xges.

    Parameters
    ----------
    data : pd.DataFrame
        The data to run the algorithm on.
    alpha : float (default=2.0)
        The alpha parameter for the BIC score.
    extended_search : bool (default=True)
        Whether to run the extended search or not. If False, the algorithm will only run XGES-0.
    ground_truth : nx.DiGraph (default=None)
        The ground truth graph, if available.
    """

    if DEBUG:
        tmp_path = os.path.join("./", "xges_" + str(uuid.uuid4()))
    else:
        tmp_path = os.path.join(gettempdir(), "xges_" + str(uuid.uuid4()))

    os.makedirs(tmp_path, exist_ok=True)

    in_file = os.path.join(tmp_path, "data.npy")
    out_file = os.path.join(tmp_path, "result.csv")
    stats_file = os.path.join(tmp_path, "stats.csv")
    truth_file = os.path.join(tmp_path, "truth.csv")
    if isinstance(data, pd.DataFrame):
        data = data.values
    data = np.ascontiguousarray(data)
    data = data.astype(np.float64)
    np.save(in_file, data)

    # 2) run xges by calling the c++ binary:
    # xges --input <data input> -a <alpha> --output <output file> --stats <stats file> -v0 [-0] [-g <ground truth>]
    # -0: run only XGES-0 (no extended search)
    xges_path = "../cmake-build-release/src/xges"
    cmd = [
        xges_path,
        "--input",
        in_file,
        "-a",
        str(alpha),
        "--output",
        out_file,
        "--stats",
        stats_file,
    ]
    if ground_truth is not None:
        df = pd.DataFrame(nx.to_numpy_array(ground_truth, nodelist=list(range(data.shape[1]))))
        df.to_csv(truth_file, index=False)
        cmd += ["-g", truth_file]
    if not extended_search:
        cmd += ["-0"]
    cmd += [f"-v{verbose}"]  # verbose level

    if DEBUG:
        print(" ".join(cmd))
        exit()

    if baseline:
        cmd += ["-b", baseline]

    error_code = subprocess.run(cmd, timeout=60 * 60 * 24)
    if error_code.returncode != 0:
        print("error running xges")
        # show the error message, the cmd used and the tmp location
        print(error_code)
        print(" ".join(cmd))
        print(tmp_path)
        raise RuntimeError("error running xges")

    # 3) load the result from the tmp location, pd.read_csv(out_file)
    graph = nx.DiGraph(pd.read_csv(out_file, header=0).values)
    xges_stats = pd.read_csv(stats_file, index_col=0, header=None).T.to_dict("records")[0]

    # 4) remove the tmp location
    if not DEBUG:
        shutil.rmtree(tmp_path)

    return graph, xges_stats


def simulation(
    d,
    n_control=1000,
    edges_per_d=2,
    seed=0,
    normalized=True,
    max_variance=1,
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
    weights = true_adjacency * rng.uniform(1, 3, true_adjacency.shape)
    if not positive:
        weights *= rng.choice([-1, 1], true_adjacency.shape)
    if normalized:
        weights = weights / (np.abs(weights).sum(axis=0, keepdims=True) + 1e-5)

    for e in true_dag.edges:
        true_dag.edges[e]["w"] = weights[e[0], e[1]]
    scm = sempler.LGANM(weights, (-2, 2), (0, max_variance), random_state=seed)
    data = scm.sample(
        n=n_control,
        random_state=seed,
    )

    data = data.astype(np.float64)
    return data, true_dag


def run_filter_kwargs(func, all_kwargs, **kwargs):
    func_args = inspect.signature(func).parameters
    all_kwargs = all_kwargs.copy()
    all_kwargs.update(kwargs)
    return func(**{k: v for k, v in all_kwargs.items() if k in func_args})


def worker(kwargs_):
    data, true_graph = run_filter_kwargs(simulation, kwargs_)
    res = kwargs_.copy()
    res_2 = run_method_given_data(data=data, true_graph=true_graph, **kwargs_)
    res.update(res_2)
    res["n_obs"] = data.shape[0]
    print(res)
    return res


def compare_methods(n_jobs=1, **kwargs):
    import multiprocessing
    from itertools import product

    if n_jobs == 1:
        scores = [
            worker(dict(zip(kwargs.keys(), kwargs_))) for kwargs_ in product(*kwargs.values())
        ]
    else:
        with multiprocessing.Pool(n_jobs) as p:
            scores = p.map(
                worker, [dict(zip(kwargs.keys(), kwargs_)) for kwargs_ in product(*kwargs.values())]
            )
    scores = pd.DataFrame(scores)
    scores = scores[sorted(scores.columns)]
    return scores


def run_method_given_data(method_name, data, true_graph, **kwargs):
    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data)

    method_func = {
        "ges": run_ges,
        "ges+turn": run_ges,
        "fges": run_fges,
        "xges0": run_xges,
        "xges": run_xges,
        "ops": run_xges,
        "b-ges": run_xges,
        "b-ges-r": run_xges,
    }[method_name]

    if method_name == "ges":
        kwargs["phases"] = ("forward", "backward")
    elif method_name == "ges+turn":
        kwargs["phases"] = ("forward", "backward", "turning")

    if method_name == "xges0":
        kwargs["extended_search"] = False
    if method_name == "xges1":
        kwargs["extended_search"] = True
    if method_name == "ops":
        kwargs["baseline"] = "ops"
    if method_name == "b-ges":
        kwargs["baseline"] = "ges"
    if method_name == "b-ges-r":
        kwargs["baseline"] = "ges-r"

    tmp_start_time = time.time()
    graph = run_filter_kwargs(method_func, kwargs, data=data, true_graph=true_graph)
    res = dict()
    if isinstance(graph, tuple):
        graph, stats = graph
        for k, v in stats.items():
            res[k] = v
        # res["score"] = stats["score"]
    res["time"] = time.time() - tmp_start_time

    # true_graph is a DAG; graph is a CPDAG where undirected edges are represented as two directed edges

    # compute metrics
    res["ecshd"] = compute_shd_cpdag(true_graph, graph)
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


DEBUG = False


def main():
    # warmup jvm and R bridge
    _ = compare_methods(
        n_jobs=1,
        seed=list(range(1)),
        n_control=[100],
        method_name=[
            "ges",
            "ges+turn",
            "fges",
            "xges0",
            "xges",
        ],
        d=[4],
        edges_per_d=[1],
        alpha=[2],
        max_variance=[0.5],
    )

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(30)),
        n_control=[100, 500, 1000, 5000, 10_000, 50_000, 100_000, 500_000, 1_000_000],
        method_name=[
            "xges0",
            "ops",
            "xges",
            "b-ges",
            "b-ges-r",
            "ges",
            "ges+turn",
            "fges",
        ],
        d=[50],
        edges_per_d=[3],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/benchmark-fig5-v2-n.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(30)),
        n_control=[10_000],
        method_name=[
            "xges0",
            # "ops-r",
            # "ops",
            "xges",
            "b-ges",
            "b-ges-r",
            "ges",
            "ges+turn",
            "fges",
        ],
        d=[50],
        edges_per_d=[3],
        alpha=[1, 2, 3, 5, 10, 15, 20, 30],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/benchmark-fig5-v2-alpha.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(30)),
        n_control=[int(n) for n in np.logspace(2, 8, 30)],
        method_name=[
            "xges0",
            "xges",
            "b-ges",
            "b-ges-r",
        ],
        d=[15],
        edges_per_d=[2],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/benchmark-double-descent.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(30)),
        n_control=[10_000],
        method_name=[
            "xges0",
            "xges",
            "ges",
            "ges+turn",
            "ops",
            "fges",
            "b-ges",
            "b-ges-r",
        ],
        d=[25, 50],
        edges_per_d=[2, 3, 4],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/benchmark-main-fig.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(10)),
        n_control=[10_000],
        method_name=["ges"],
        d=[10, 50, 100],
        edges_per_d=[1, 2, 3, 4],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/study-ges.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(5)),
        n_control=[10_000],
        method_name=[
            "xges0",
            "xges",
            "b-ges",
            "b-ges-r",
        ],
        d=[10, 25, 50, 100, 200, 300, 500],
        edges_per_d=[2, 4],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/benchmark-rebuttal-n_bic.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(5)),
        n_control=[10_000],
        method_name=[
            "xges0",
            "xges",
            "ges",
            "ges+turn",
            "ops",
            "fges",
        ],
        d=list(np.logspace(1, 2.75, 8).astype(int)),
        edges_per_d=[1, 2, 3, 4],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/vary-d-and-rho-seeds-1-to-5.csv", index=False)

    start = time.time()
    res = compare_methods(
        n_jobs=1,
        seed=list(range(5)),
        n_control=[10_000],
        method_name=[
            "xges0",
            "xges",
            "ges",
            "ges+turn",
            "fges",
        ],
        d=[10, 50, 125, 250, 500, 750],
        edges_per_d=[2, 4],
        alpha=[2],
        n_inter=[0],
        max_variance=[0.5],
    )
    total_time = time.time() - start
    print("total time:", total_time)
    res.to_csv("results/benchmark-speed.csv", index=False)


if __name__ == "__main__":
    # warmup jvm and R bridge
    _ = compare_methods(
        n_jobs=1,
        seed=list(range(1)),
        n_control=[100],
        method_name=[
            "ges",
            "ges+turn",
            "fges",
            "xges0",
            "xges",
        ],
        d=[4],
        edges_per_d=[1],
        alpha=[2],
        max_variance=[0.5],
    )
    main()
