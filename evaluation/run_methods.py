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

    if intervention:
        data = [data]
        intervention_labels = [-1] * n_control
        for i in range(d):
            intervention_labels += [i] * n_inter
            data.append(scm.sample(n=n_inter, random_state=i, shift_interventions={i: (-0.1, 0)}))
        data = np.concatenate(data, axis=0)
        intervention_labels = np.array(intervention_labels)

    # normalize the data? only on the control data?

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


def run_xges(data, alpha, **kwargs):
    # 1) save the data to a tmp location, generate a unique tmp file name (e.g. using uuid)

    tmp_path = os.path.join(gettempdir(), "xges_" + str(uuid.uuid4()))
    os.makedirs(tmp_path, exist_ok=True)
    in_file = os.path.join(tmp_path, "data.npy")
    out_file = os.path.join(tmp_path, "result.csv")
    stats_file = os.path.join(tmp_path, "stats.csv")
    np.save(in_file, data)

    # 2) run xges by calling the c++ binary: xges -f <tmp location> -a <alpha>
    xges_path = "../cmake-build-release/ges_cpp"
    cmd = [xges_path, "-i", in_file, "-a", str(alpha), "-o", out_file, "--stats", stats_file]
    start_time = time.time()
    subprocess.run(cmd)
    xges_time = time.time() - start_time

    # 3) load the result from the tmp location, pd.read_csv(out_file)
    graph = nx.DiGraph(pd.read_csv(out_file, header=0).values)
    xges_stats = pd.read_csv(stats_file, index_col=0, header=None).T.to_dict("records")[0]
    print(xges_stats)
    xges_stats["time_python"] = xges_time

    # 4) remove the tmp location
    shutil.rmtree(tmp_path)

    return graph, xges_stats


from pdag import compute_ecshd_from_graph, PDAG


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
                n_inter=None,
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


if __name__ == "__main__":
    benchmark_xges_alpha()
    # TODO: just call run_xges
