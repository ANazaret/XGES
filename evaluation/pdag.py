import copy
from collections import defaultdict

import networkx as nx


class PDAG:
    """
    Partial Directed Acyclic Graph (PDAG).

    A directed graph that may contain undirected edges.
    """

    def __init__(self, nodes):
        self.nodes = set(nodes)
        self.children = defaultdict(set)
        self.parents = defaultdict(set)
        self.neighbors = defaultdict(set)
        self.adjacent = defaultdict(set)

        self.number_of_undirected_edges = 0
        self.number_of_directed_edges = 0

    @property
    def number_of_edges(self):
        return self.number_of_directed_edges + self.number_of_undirected_edges

    def __eq__(self, other):
        same_directed_edges = all([self.children[node] == other.children[node] for node in self.nodes])
        same_undirected_edges = all([self.neighbors[node] == other.neighbors[node] for node in self.nodes])
        return same_directed_edges and same_undirected_edges

    def get_directed_edges(self):
        edges = set()
        for node in self.nodes:
            for child in self.children[node]:
                edges.add((node, child))
        return edges

    def get_undirected_edges(self):
        edges = set()
        for node in self.nodes:
            for neighbor in self.neighbors[node]:
                if node < neighbor:
                    edges.add((node, neighbor))
        return edges

    def get_parents(self, node) -> set:
        return self.parents[node]

    def get_children(self, node) -> set:
        return self.children[node]

    def get_neighbors(self, node) -> set:
        return self.neighbors[node]

    def get_adjacent(self, node) -> set:
        res = self.adjacent[node]
        return res

    def has_directed_edge(self, x, y):
        return y in self.children[x]

    def has_undirected_edge(self, x, y):
        return y in self.neighbors[x]

    def remove_directed_edge(self, x, y):
        self.children[x].remove(y)
        self.parents[y].remove(x)
        self.adjacent[x].remove(y)
        self.adjacent[y].remove(x)
        self.number_of_directed_edges -= 1

    def remove_undirected_edge(self, x, y):
        self.neighbors[x].remove(y)
        self.neighbors[y].remove(x)
        self.adjacent[x].remove(y)
        self.adjacent[y].remove(x)
        self.number_of_undirected_edges -= 1

    def add_directed_edge(self, x, y):
        self.children[x].add(y)
        self.parents[y].add(x)
        self.adjacent[x].add(y)
        self.adjacent[y].add(x)
        self.number_of_directed_edges += 1

    def add_undirected_edge(self, x, y):
        self.neighbors[x].add(y)
        self.neighbors[y].add(x)
        self.adjacent[x].add(y)
        self.adjacent[y].add(x)
        self.number_of_undirected_edges += 1

    def get_edge_status(self, x, y):
        if self.has_directed_edge(x, y):
            return "-->"
        if self.has_directed_edge(y, x):
            return "<--"
        if self.has_undirected_edge(x, y):
            return "---"
        return "xxx"

    def __repr__(self):
        edges = [f"{u} -> {v}" for u in sorted(self.children) for v in sorted(self.children[u])] + [
            f"{u} -- {v}" for u in sorted(self.neighbors) for v in sorted(self.neighbors[u]) if u < v
        ]
        return f"PDAG({', '.join(edges)})"

    def remove_node(self, node):
        self.nodes.remove(node)
        for u in list(self.neighbors[node]):
            self.remove_undirected_edge(node, u)
        for u in list(self.children[node]):
            self.remove_directed_edge(node, u)
        for u in list(self.parents[node]):
            self.remove_directed_edge(u, node)

    def to_dag(self):
        return pdag_to_dag(self)

    def to_cpdag(self):
        dag = self.to_dag()
        return dag_to_cpdag(dag)

    def shd_against_dag(self, dag: nx.DiGraph):
        """
        Returns the structural hamming distance between the PDAG and a DAG, where if an edge is undirected in the PDAG
        but directed in the DAG, it counts as 0.
        """
        shd = 0
        for x in self.nodes:
            for y in self.nodes:
                if x >= y:
                    continue
                if self.get_edge_status(x, y) == "---":
                    if get_dag_edge_status(dag, x, y) == "xxx":
                        shd += 1
                elif self.get_edge_status(x, y) != get_dag_edge_status(dag, x, y):
                    shd += 1
        return shd

    @staticmethod
    def shd_cpdag(pdag1, pdag2):
        score = 0
        # for each pair of nodes (x,y)
        # we either have: x  y; x -- y; x -> y; x <- y.
        for x in pdag1.nodes:
            for y in pdag1.nodes:
                if x >= y:
                    continue
                if pdag1.get_edge_status(x, y) != pdag2.get_edge_status(x, y):
                    score += 1
        return score

    @staticmethod
    def from_digraph(dag):
        pdag = PDAG(dag.nodes)
        for x, y in dag.edges:
            if dag.has_edge(y, x):
                if not pdag.has_undirected_edge(x, y):
                    pdag.add_undirected_edge(x, y)
            else:
                pdag.add_directed_edge(x, y)
        return pdag


def get_dag_edge_status(dag: nx.DiGraph, x, y):
    if dag.has_edge(x, y) and dag.has_edge(y, x):
        raise ValueError("Edge is undirected in the DAG. Not a valid DAG.")
    if dag.has_edge(x, y):
        return "-->"
    if dag.has_edge(y, x):
        return "<--"
    return "xxx"


def pdag_to_dag(pdag: PDAG):
    """
    Returns a consistent DAG from the PDAG.
    Algorithm [1] as described in [2].

    References
    ----------
    [1] https://ftp.cs.ucla.edu/pub/stat_ser/r185-dor-tarsi.pdf
    [2] https://www.jmlr.org/papers/volume3/chickering02b/chickering02b.pdf
    """

    output_dag = nx.DiGraph()
    output_dag.add_nodes_from(pdag.nodes)
    pdag = copy.deepcopy(pdag)

    # while the pdag is not empty
    while pdag.nodes:
        # find a node x that:
        # 1. is a sink (no directed outgoing edges)
        # 2. each neighbor y of x is adjacent to all other adjacent nodes of x
        # Note: as pointed out by [github ges], this is not the same as saying it is a clique [2] (which is not correct)
        found_valid_x = False
        for x in pdag.nodes:
            # condition 1: sink
            if len(pdag.get_children(x)):
                continue

            # condition 2: each neighbor y of x is adjacent to all other adjacent nodes of x
            ne_x = pdag.get_neighbors(x)
            ad_x = pdag.get_adjacent(x)
            is_valid = True
            for y in ne_x:
                if not ad_x.issubset(pdag.get_adjacent(y) | {y}):
                    is_valid = False
                    break
            if not is_valid:
                continue

            # node x is valid
            found_valid_x = True
            # each undirected edge between x and y is directed as y → x in the output DAG
            for y in ne_x:
                output_dag.add_edge(y, x)
            for y in pdag.get_children(x):  # seems useless
                output_dag.add_edge(x, y)
            for y in pdag.get_parents(x):
                output_dag.add_edge(y, x)
            # remove x from the PDAG
            pdag.remove_node(x)
            break

        if not found_valid_x:
            raise ValueError("No valid node x found, PDAG is not consistent.")

    return output_dag


def _order_edges(dag: nx.DiGraph) -> list:
    """
    Returns a total ordering of the edges in a DAG.

    Algorithm Order-Edges from [1, Figure 13].
    """
    node_ordering = list(nx.topological_sort(dag))
    node_to_order = {node: i for i, node in enumerate(node_ordering)}
    ordered_edges = []
    for node in node_ordering:
        predecessors = list(dag.predecessors(node))
        # sort predecessors by order
        predecessors.sort(key=lambda x: node_to_order[x], reverse=True)
        ordered_edges.extend([(p, node) for p in predecessors])
    return ordered_edges


def dag_to_cpdag(dag) -> PDAG:
    """
    Returns an essential graph from a DAG.
    Algorithm Label-Edges from [1, Figure 14].

    References
    ----------
    [1] https://www.jmlr.org/papers/volume3/chickering02b/chickering02b.pdf
    """
    ordered_edges = _order_edges(dag)
    pdag = PDAG(dag.nodes)
    for x, y in ordered_edges:
        if pdag.has_undirected_edge(x, y) or pdag.has_directed_edge(x, y):
            # edge is already labeled
            continue
        break_loop = False
        for w in pdag.get_parents(x):
            if not dag.has_edge(w, y):
                # label x -> y and all edges incident into y as compelled
                # pdag.directed_graph.add_edge(x, y)
                for z in dag.predecessors(y):
                    pdag.add_directed_edge(z, y)
                break_loop = True
                break
            else:
                # label w -> y as compelled
                pdag.add_directed_edge(w, y)
        if break_loop:
            continue
        # if there exists an edge z -> y such that z != x and z is not a parent of x
        if set(dag.predecessors(y)) - {x} - set(dag.predecessors(x)):
            # label x -> y and all edges incident into y as compelled
            for z in dag.predecessors(y):
                if pdag.has_undirected_edge(z, y):
                    continue
                pdag.add_directed_edge(z, y)
        else:
            # label x -> y and all edges incident into y as reversible
            for z in dag.predecessors(y):
                if pdag.has_directed_edge(z, y):
                    continue
                pdag.add_undirected_edge(z, y)

    return pdag


def compute_shd_cpdag(graph1: nx.DiGraph, graph2: nx.DiGraph) -> int:
    cpdag1 = PDAG.from_digraph(graph1).to_cpdag()
    cpdag2 = PDAG.from_digraph(graph2).to_cpdag()

    return PDAG.shd_cpdag(cpdag1, cpdag2)
