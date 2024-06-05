import time
from collections import defaultdict, deque

import numpy as np


class PDAG:
    """
    Represents a partially directed acyclic graph (PDAG).

    A PDAG is a graph that contains both directed and undirected edges.

    Parameters
    ----------
    num_nodes : int
        The number of nodes in the graph.

    """

    # TODO replace defaultdict with list of sets
    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self.children = defaultdict(set)
        self.parents = defaultdict(set)
        self.neighbors = defaultdict(set)
        self.adjacent = defaultdict(set)
        self.adjacent_reachable = defaultdict(set)

        self.number_of_undirected_edges = 0
        self.number_of_directed_edges = 0

        self.nodes = set(range(self.num_nodes))

        self.forbidden_insert_parents = defaultdict(set)

        self.block_semi_directed_path_visited = np.zeros(self.num_nodes, dtype=int)
        self.block_semi_directed_path_blocked = np.zeros(self.num_nodes, dtype=int)
        self.block_semi_directed_path_parent = np.zeros(self.num_nodes, dtype=int)
        self.block_semi_directed_path_queue = []

        self.statistics = defaultdict(int)

    def get_nodes(self):
        return self.nodes

    def get_number_of_edges(self):
        return self.number_of_directed_edges + self.number_of_undirected_edges

    def get_parents(self, node):
        return self.parents[node]

    def get_children(self, node):
        return self.children[node]

    def get_neighbors(self, node):
        return self.neighbors[node]

    def get_adjacent(self, node):
        return self.adjacent[node]

    def get_adjacent_reachable(self, node):
        return self.adjacent_reachable[node]

    def get_neighbors_adjacent(self, node_y, node_x):
        return self.neighbors[node_y] & self.adjacent[node_x]

    def get_neighbors_not_adjacent(self, node_y, node_x):
        return self.neighbors[node_y] - self.adjacent[node_x]

    def has_directed_edge(self, x, y):
        return y in self.children[x]

    def has_undirected_edge(self, x, y):
        return y in self.neighbors[x]

    def remove_directed_edge(self, x, y):
        self.children[x].remove(y)
        self.parents[y].remove(x)
        self.adjacent[x].remove(y)
        self.adjacent[y].remove(x)
        self.adjacent_reachable[x].remove(y)
        self.number_of_directed_edges -= 1

    def remove_undirected_edge(self, x, y):
        self.neighbors[x].remove(y)
        self.neighbors[y].remove(x)
        self.adjacent[x].remove(y)
        self.adjacent[y].remove(x)
        self.adjacent_reachable[x].remove(y)
        self.adjacent_reachable[y].remove(x)
        self.number_of_undirected_edges -= 1

    def add_directed_edge(self, x, y):
        self.children[x].add(y)
        self.parents[y].add(x)
        self.adjacent[x].add(y)
        self.adjacent[y].add(x)
        self.adjacent_reachable[x].add(y)
        self.number_of_directed_edges += 1

    def add_undirected_edge(self, x, y):
        self.neighbors[x].add(y)
        self.neighbors[y].add(x)
        self.adjacent[x].add(y)
        self.adjacent[y].add(x)
        self.adjacent_reachable[x].add(y)
        self.adjacent_reachable[y].add(x)
        self.number_of_undirected_edges += 1

    def is_clique(self, nodes: set):
        for node in nodes:
            if not nodes.issubset(self.neighbors[node]):
                return False
        return True

    def is_insert_valid(self, insert, unblocked_paths_map, reverse=False):
        start_time = time.time()
        self.statistics["is_insert_valid-#calls"] += 1
        x = insert.x
        y = insert.y
        T = insert.T

        # 0. check if edge is forbidden
        if self.is_insert_forbidden(x, y):
            self.statistics["is_insert_valid-false_0-#"] += 1
            self.statistics["is_insert_valid-false-time"] += time.time() - start_time
            return False

        adjacent_x = self.adjacent[x]
        if not reverse:
            # 1. x and y are not adjacent
            if y in adjacent_x:
                self.statistics["is_insert_valid-false_1a-#"] += 1
                self.statistics["is_insert_valid-false-time"] += time.time() - start_time
                return False
        else:
            # 1. x ← y
            parents_x = self.parents[x]
            if y not in parents_x:
                self.statistics["is_insert_valid-false_1b-#"] += 1
                self.statistics["is_insert_valid-false-time"] += time.time() - start_time
                return False

        # 2. T ⊆ Ne(y) \ Ad(x)
        # <=> T ⊆ Ne(y) and T does not intersect Ad(x)
        neighbors_y = self.neighbors[y]
        if not (T.issubset(neighbors_y) and T.isdisjoint(adjacent_x)):
            self.statistics["is_insert_valid-false_2-#"] += 1
            self.statistics["is_insert_valid-false-time"] += time.time() - start_time
            return False

        # 5. E (insert.effective_parents) == [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
        ne_y_ad_x_T = neighbors_y.intersection(adjacent_x).union(T)
        if insert.effective_parents != ne_y_ad_x_T.union(self.parents[y]):
            self.statistics["is_insert_valid-false_5-#"] += 1
            self.statistics["is_insert_valid-false-time"] += time.time() - start_time
            return False

        # 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
        if not self.is_clique(ne_y_ad_x_T):
            self.statistics["is_insert_valid-false_3-#"] += 1
            self.statistics["is_insert_valid-false-time"] += time.time() - start_time
            return False

        # 4. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x
        ignore_direct_edge = reverse
        if reverse:
            # for reverse: ne_y_ad_x_T is actually [Ne(y) ∩ Ad(x)] ∪ T ∪ Ne(x)
            ne_y_ad_x_T = ne_y_ad_x_T.union(self.neighbors[x])
        if not self.is_blocking_semi_directed_paths(
                y, x, ne_y_ad_x_T, unblocked_paths_map, ignore_direct_edge
        ):
            self.statistics["is_insert_valid-false_4-#"] += 1
            self.statistics["is_insert_valid-false-time"] += time.time() - start_time
            return False

        self.statistics["is_insert_valid-true-#"] += 1
        self.statistics["is_insert_valid-true-time"] += time.time() - start_time
        return True

    def is_reverse_valid(self, reverse, unblocked_paths_map):
        # is Pa(x) unchanged
        if self.parents[reverse.insert.x] != reverse.parents_x:
            return False
        return self.is_insert_valid(reverse.insert, unblocked_paths_map, True)

    def is_delete_valid(self, delet):
        # 1. x and y are neighbors or x is a parent of y [aka y is adjacent_reachable from x]
        x = delet.x
        y = delet.y
        if y not in self.adjacent_reachable[x]:
            return False

        # 2. C is a subset of Ne(y) ∩ Ad(x)
        # <=> C ⊆ Ne(y) and C ⊆ Ad(x)
        neighbors_y = self.neighbors[y]
        adjacent_x = self.adjacent[x]
        if not (delet.C.issubset(neighbors_y) and delet.C.issubset(adjacent_x)):
            return False

        # 3. E (delet.effective_parents) = C ∪ Pa(y) ∪ {x}
        if delet.effective_parents != delet.C.union(self.parents[y]).union({x}):
            return False

        # 4. C is a clique
        if not self.is_clique(delet.C):
            return False

        return True

    def is_insert_forbidden(self, x, y):
        return x in self.forbidden_insert_parents[y]

    def is_blocking_semi_directed_paths(self, y, x, blocked_nodes, unblocked_paths_map, ignore_direct_edge):
        if y == x:
            return False
        self.statistics["block_semi_directed_paths-#calls"] += 1
        start_time = time.time()

        # Set block_semi_directed_path_visited to 0
        self.block_semi_directed_path_visited.fill(0)
        visited = self.block_semi_directed_path_visited
        self.block_semi_directed_path_blocked.fill(0)
        blocked = self.block_semi_directed_path_blocked
        for n in blocked_nodes:
            blocked[n] = 1

        # BFS search from y to x, using adjacent_reachable edges, avoiding blocked nodes
        visited[y] = 1

        queue = deque()
        queue.append(y)

        while queue:
            node = queue.popleft()
            reachable = self.get_adjacent_reachable(node)

            for n in reachable:
                if visited[n] or blocked[n] or (node == y and n == x and ignore_direct_edge):
                    continue
                self.block_semi_directed_path_parent[n] = node
                if n == x:
                    self.statistics["block_semi_directed_paths-false-#"] += 1
                    self.statistics["block_semi_directed_paths-false-time"] += time() - start_time
                    # retrieve the path
                    current = x
                    while current != y:
                        parent = self.block_semi_directed_path_parent[current]
                        unblocked_paths_map[(parent, current)].add((x, y))
                        current = parent
                    return False
                queue.append(n)
                visited[n] = 1

        self.statistics["block_semi_directed_paths-true-#"] += 1
        self.statistics["block_semi_directed_paths-true-time"] += time.time() - start_time
        return True
