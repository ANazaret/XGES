import numpy as np
import pytest

from xges import XGES


def generate_independent_data(n, d):
    X = np.random.randn(n, d)
    return X


@pytest.fixture()
def generate_v_structure_data(n=200):
    X = np.random.randn(n, 3)
    X[:, 1] = X[:, 1] * 0.5 + X[:, 0] + X[:, 2]
    return X


@pytest.fixture()
def generate_chain_data(n=200):
    X = np.random.randn(n, 3)
    X[:, 1] = X[:, 0] + X[:, 1] * 0.5
    X[:, 2] = X[:, 1] + X[:, 2] * 0.5
    return X


def test_trivial_data():
    data = generate_independent_data(200, 3)

    xges = XGES()
    xges.fit(data, use_fast_numba=False)
    pdag = xges.get_pdag()
    assert pdag.get_number_of_edges() == 0


@pytest.mark.parametrize(
    "data, directed_edges, undirected_edges",
    [
        ("generate_v_structure_data", [(0, 1), (2, 1)], []),
        ("generate_chain_data", [], [(0, 1), (1, 2)]),
    ],
)
def test_three_nodes(data, directed_edges, undirected_edges, request):
    data = request.getfixturevalue(data)
    xges = XGES()
    xges.fit(data, use_fast_numba=False)
    pdag = xges.get_pdag()

    assert pdag.number_of_directed_edges == len(directed_edges)
    assert pdag.number_of_undirected_edges == len(undirected_edges)
    for edge in directed_edges:
        assert pdag.has_directed_edge(*edge)

    for edge in undirected_edges:
        assert pdag.has_undirected_edge(*edge)


def test_numba():
    data = generate_independent_data(200, 3)

    xges = XGES()
    xges.fit(data, use_fast_numba=True)
    pdag = xges.get_pdag()
    assert pdag.get_number_of_edges() == 0
