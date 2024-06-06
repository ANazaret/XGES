# Extremely Greedy Equivalence Search

Extremely Greedy Equivalence Search (XGES) is an efficient algorithm
for learning the structure of a causal graph from observational data.
It improves upon the Greedy Equivalence Search (GES) algorithm with a
more accurate search strategy and a more efficient implementation.

The algorithm is described in the paper "Extremely Greedy Equivalence Search" (see [citation below](#Citation))

In this repo, we provide:

- a pure python implementation of the algorithm in `xges/` available with `pip install xges`.
    - it can optionally use `numba` for faster execution.
    - examples of usage are available in `examples/simple.py`.
- a pure c++ implementation of the algorithm in `src-cpp/`, which is at least ~10x faster
  than the python implementation (with `numba`).
- code to reproduce the experiments in the paper in `evaluation/`.
    - an ad-hoc python wrapper calling the cpp executable is available in `evaluation/benchmarks.py`.
    - a notebook to generate the figures in the paper in `evaluation/paper.ipynb`.

## Using the python package

The python package can be installed with pip:

```bash
pip install xges
```

The package can be used as follows:

```python
from xges import XGES

data = ...
xges = XGES()
pdag = xges.fit(data)  # PDAG object representing the Markov equivalence class (MEC)

# PDAG object with only directed edges, representing an arbitrary DAG in the MEC
a_dag = pdag.get_dag_extension()

# networkx DiGraph object with two edges for undirected PDAG edges
networkx_pdag = pdag.to_networkx()
networkx_dag = a_dag.to_networkx()

adjacency_pdag = pdag.to_adjacency_matrix()
adjacency_dag = a_dag.to_adjacency_matrix()



```

## Reproducing the experiments

The experiments can be reproduced by running the `evaluation/benchmarks.py`
script (after compiling the c++ code in `src-cpp`).
The figures are generated in the notebook `evaluation/paper.ipynb`.

## Building the c++ code

Use the CMakeLists.txt file to build the code.

The code can be run with the following command:

```bash
xges --input data.npy --output out.csv --stats stats.csv -v1
```

The input file should be a numpy file with the data matrix. The output file
will contain the CPDAG. The stats file will contain some statistics collected
during the execution of the algorithm.
`-v1` is the verbosity level. It can be set to 0, 1, or 2.

More options can be found by running `xges --help`.

## Citation

If you use this code, please cite the following paper:

```
@inproceedings{nazaret2021extremely,
  title={Extremely Greedy Equivalence Search},
  author={Nazaret, Achille and Blei, David},
  booktitle={Uncertainty in Artificial Intelligence},
  year={2024}
}
```