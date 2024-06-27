import numpy as np

from xges import XGES

# Example 1: independent data
n = 200
d = 3
data = np.random.normal(size=(n, d))
model = XGES()
model.fit(data)
print("Independent:", model.get_pdag())

# Example 2: v-structure
n = 200
data = np.random.normal(size=(n, 3))
data[:, 1] = data[:, 1] * 0.5 + data[:, 0] + data[:, 2]
model = XGES()
model.fit(data)
print("V-structure:", model.get_pdag())

# Example 3: chain
n = 200
data = np.random.normal(size=(n, 3))
data[:, 1] = data[:, 0] + data[:, 1] * 0.5
data[:, 2] = data[:, 1] + data[:, 2] * 0.5
model = XGES()
model.fit(data)
print("Chain:", model.get_pdag())

# Get one possible DAG
print("Chain, one possible DAG:", model.get_a_dag())

# Get the DAG as a networkx.DiGraph if networkx is installed
try:
    dag_nx = model.get_a_dag().to_networkx()

    print("Chain, one possible DAG as a networkx.DiGraph:", dag_nx)
except ImportError:
    dag_nx = None
    print("networkx not installed.")

# Get the adjacency matrix of the PDAG
print("Chain, adjacency matrix:", model.get_pdag().to_adjacency_matrix())
