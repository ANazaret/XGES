# Extremely Greedy Equivalence Search

This is the code for the paper "Extremely Greedy Equivalence Search".

We recommend checking any updated code at
the [official repository](https://github.com/ANazaret/XGES).

## Reproducing the experiments

The experiments can be reproduced by running the `evaluation/benchmarks.py`
script.
The figures are generated in the notebook `evaluation/paper.ipynb`.

## Building the code

Use the CMakeLists.txt file to build the code.

## Running the code

We recommend checking the simple Python wrapper `evaluation/benchmarks.py` to
see how to call the `xges` executable.

The code can be run with the following command:

```bash
xges --input data.npy --output out.csv --stats stats.csv -v1
```

The input file should be a numpy file with the data matrix. The output file
will contain the CPDAG. The stats file will contain some statistics collected
during the execution of the algorithm.
`-v1` is the verbosity level. It can be set to 0, 1, or 2.

More options can be found by running `xges --help`.

## Citing

If you use this code, please cite the following paper:

```
@inproceedings{nazaret2021extremely,
  title={Extremely Greedy Equivalence Search},
  author={Nazaret, Achille and Blei, David},
  booktitle={Uncertainty in Artificial Intelligence},
  year={2024}
}
```