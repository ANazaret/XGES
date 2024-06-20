import copy
import time
from collections import defaultdict, deque

import numpy as np

from xges.edge_queue_set import EdgeQueueSet, EdgeType


class PDAG:
    """
    Represents a partially directed acyclic graph (PDAG).

    A PDAG is a graph that contains both directed and undirected edges.

    Parameters
    ----------
    num_nodes : int
        The number of nodes in the graph.

    """

    def __init__(self, num_nodes):
        self.num_nodes = num_nodes
        self.children = [set() for _ in range(self.num_nodes)]
        self.parents = [set() for _ in range(self.num_nodes)]
        self.neighbors = [set() for _ in range(self.num_nodes)]
        self.adjacent = [set() for _ in range(self.num_nodes)]
        self.adjacent_reachable = [set() for _ in range(self.num_nodes)]

        self.number_of_undirected_edges = 0
        self.number_of_directed_edges = 0

        self.nodes = set(range(self.num_nodes))

        self.forbidden_insert_parents = [set() for _ in range(self.num_nodes)]

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

    def is_empty(self):
        # To change if we allow known edges at the beginning
        return self.number_of_directed_edges == 0 and self.number_of_undirected_edges == 0

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

    def add_forbidden_insert(self, x, y):
        self.forbidden_insert_parents[y].add(x)

    def apply_edge_modification(self, edge_modification, undo=False):
        old_type = edge_modification.old_type
        new_type = edge_modification.new_type
        if undo:
            old_type, new_type = new_type, old_type
        if old_type == EdgeType.UNDIRECTED:
            self.remove_undirected_edge(edge_modification.x, edge_modification.y)
        elif old_type == EdgeType.DIRECTED_TO_X:
            self.remove_directed_edge(edge_modification.y, edge_modification.x)
        elif old_type == EdgeType.DIRECTED_TO_Y:
            self.remove_directed_edge(edge_modification.x, edge_modification.y)
        if new_type == EdgeType.UNDIRECTED:
            self.add_undirected_edge(edge_modification.x, edge_modification.y)
        elif new_type == EdgeType.DIRECTED_TO_X:
            self.add_directed_edge(edge_modification.y, edge_modification.x)
        elif new_type == EdgeType.DIRECTED_TO_Y:
            self.add_directed_edge(edge_modification.x, edge_modification.y)

    def is_clique(self, nodes: set):
        for node in nodes:
            adjacent = self.get_adjacent(node).copy()
            adjacent.add(node)
            if not nodes.issubset(adjacent):
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

    def is_delete_valid(self, delete):
        # 1. x and y are neighbors or x is a parent of y [aka y is adjacent_reachable from x]
        x = delete.x
        y = delete.y
        if y not in self.adjacent_reachable[x]:
            return False

        # 2. C is a subset of Ne(y) ∩ Ad(x)
        # <=> C ⊆ Ne(y) and C ⊆ Ad(x)
        neighbors_y = self.neighbors[y]
        adjacent_x = self.adjacent[x]
        if not (delete.C.issubset(neighbors_y) and delete.C.issubset(adjacent_x)):
            return False

        # 3. E (delete.effective_parents) = C ∪ Pa(y) ∪ {x}
        if delete.effective_parents != delete.C.union(self.parents[y]).union({x}):
            return False

        # 4. C is a clique
        if not self.is_clique(delete.C):
            return False

        return True

    def is_insert_forbidden(self, x, y):
        return x in self.forbidden_insert_parents[y]

    def is_blocking_semi_directed_paths(
        self, y, x, blocked_nodes, unblocked_paths_map, ignore_direct_edge
    ):
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
                    self.statistics["block_semi_directed_paths-false-time"] += (
                        time.time() - start_time
                    )
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

    def apply_insert(self, insert, edge_modifications_map):
        start_time = time.time()
        x = insert.x
        y = insert.y
        T = insert.T

        edges_to_check = EdgeQueueSet()

        # a. insert the directed edge x → y
        self.add_directed_edge(x, y)
        edge_modifications_map.update_edge_directed(x, y, EdgeType.NONE)
        # b. for each t ∈ T: orient the (previously undirected) edge between t and y as t → y
        for t in T:
            self.remove_undirected_edge(t, y)
            self.add_directed_edge(t, y)
            edge_modifications_map.update_edge_directed(t, y, EdgeType.UNDIRECTED)

        # check if the orientation of the surrounding edges should be updated
        edges_to_check.push_directed(x, y)
        self.add_adjacent_edges(x, y, edges_to_check)
        # edges t → y are part of a v-structure with x, so we don't need to check them
        for t in T:
            self.add_adjacent_edges(t, y, edges_to_check)

        self.complete_cpdag_efficient(edges_to_check, edge_modifications_map)
        self.statistics["apply_insert-time"] += time.time() - start_time

    def apply_reverse(self, reverse, edge_modifications_map):
        x = reverse.insert.x
        y = reverse.insert.y
        # 1. remove the directed edge y → x
        self.remove_directed_edge(y, x)
        edge_modifications_map.update_edge_none(x, y, EdgeType.DIRECTED_TO_X)
        # 2. apply the insert
        self.apply_insert(reverse.insert, edge_modifications_map)

    def apply_delete(self, delet, edge_modifications_map):
        start_time = time.time()
        if self.has_directed_edge(delet.x, delet.y):
            # 1. remove the directed edge x → y
            self.remove_directed_edge(delet.x, delet.y)
            edge_modifications_map.update_edge_none(delet.x, delet.y, EdgeType.DIRECTED_TO_Y)
        else:
            # 1. remove the undirected edge x - y
            self.remove_undirected_edge(delet.x, delet.y)
            edge_modifications_map.update_edge_none(delet.x, delet.y, EdgeType.UNDIRECTED)

        # H = Ne(y) ∩ Ad(x) \ C
        H = self.get_neighbors(delet.y).intersection(self.get_adjacent(delet.x)).difference(delet.C)

        # 2. for each h ∈ H:
        #   - orient the (previously undirected) edges between h and y as y → h
        #   - orient any (previously undirected) edges between x and h as x → h
        for h in H:
            self.remove_undirected_edge(delet.y, h)
            self.add_directed_edge(delet.y, h)
            edge_modifications_map.update_edge_directed(delet.y, h, EdgeType.UNDIRECTED)

            if self.has_undirected_edge(delet.x, h):
                self.remove_undirected_edge(delet.x, h)
                self.add_directed_edge(delet.x, h)
                edge_modifications_map.update_edge_directed(delet.x, h, EdgeType.UNDIRECTED)

        edges_to_check = EdgeQueueSet()
        self.add_adjacent_edges(delet.x, delet.y, edges_to_check)
        for h in H:
            self.add_adjacent_edges(h, delet.y, edges_to_check)
            self.add_adjacent_edges(delet.x, h, edges_to_check)

        self.complete_cpdag_efficient(edges_to_check, edge_modifications_map)
        self.statistics["apply_delete-time"] += time.time() - start_time

    def add_adjacent_edges(self, x, y, edge_queue_set):
        for z in self.children[y]:
            if x != z:
                edge_queue_set.push_directed(y, z)
        for z in self.parents[y]:
            if x != z:
                edge_queue_set.push_directed(z, y)
        for z in self.neighbors[y]:
            if x != z:
                edge_queue_set.push_undirected(y, z)
        for z in self.children[x]:
            if y != z:
                edge_queue_set.push_directed(x, z)
        for z in self.parents[x]:
            if y != z:
                edge_queue_set.push_directed(z, x)
        for z in self.neighbors[x]:
            if y != z:
                edge_queue_set.push_undirected(x, z)

    def is_oriented_by_meek_rule_1(self, x, y):
        for z in self.parents[x]:
            if y not in self.adjacent[z]:
                return True
        return False

    def is_oriented_by_meek_rule_2(self, x, y):
        for z in self.children[x]:
            if y in self.children[z]:
                return True
        return False

    def is_oriented_by_meek_rule_3(self, x, y):
        candidates_z_w = self.neighbors[x].intersection(self.parents[y])
        for candidate_z in candidates_z_w:
            for candidate_w in candidates_z_w:
                if candidate_z >= candidate_w:
                    continue
                if candidate_w not in self.adjacent[candidate_z]:
                    return True
        return False

    def is_oriented_by_meek_rule_4(self, x, y):
        candidates_w = self.neighbors[x].intersection(self.neighbors[y])
        for candidate_w in candidates_w:
            for candidate_z in self.children[candidate_w]:
                if y in self.children[candidate_z] and x not in self.adjacent[candidate_z]:
                    return True
        return False

    def is_part_of_v_structure(self, x, y):
        if x not in self.parents[y]:
            return False
        for z in self.parents[y]:
            if z != x and x not in self.adjacent[z]:
                return True
        return False

    def complete_cpdag_efficient(self, edges_to_check, edge_modifications_map):
        while not edges_to_check.empty():
            edge = edges_to_check.pop()
            if edge.is_directed():
                x = edge.get_source()
                y = edge.get_target()
                if not (
                    self.is_part_of_v_structure(x, y)
                    or self.is_oriented_by_meek_rule_1(x, y)
                    or self.is_oriented_by_meek_rule_2(x, y)
                    or self.is_oriented_by_meek_rule_3(x, y)
                ):
                    self.remove_directed_edge(x, y)
                    self.add_undirected_edge(x, y)
                    edge_modifications_map.update_edge_undirected(x, y, EdgeType.DIRECTED_TO_Y)
                    self.add_adjacent_edges(x, y, edges_to_check)
            else:
                x = edge.get_x()
                y = edge.get_y()
                if (
                    self.is_oriented_by_meek_rule_1(x, y)
                    or self.is_oriented_by_meek_rule_2(x, y)
                    or self.is_oriented_by_meek_rule_3(x, y)
                ):
                    pass
                elif (
                    self.is_oriented_by_meek_rule_1(y, x)
                    or self.is_oriented_by_meek_rule_2(y, x)
                    or self.is_oriented_by_meek_rule_3(y, x)
                ):
                    x, y = y, x
                else:
                    continue
                self.remove_undirected_edge(x, y)
                self.add_directed_edge(x, y)
                edge_modifications_map.update_edge_directed(x, y, EdgeType.UNDIRECTED)
                self.add_adjacent_edges(x, y, edges_to_check)
        assert self.check_is_cpdag(), "This should not happen! Report a bug please."

    def check_is_cpdag(self):
        for x in self.nodes:
            for y in self.children[x]:
                if not (
                    self.is_part_of_v_structure(x, y)
                    or self.is_oriented_by_meek_rule_1(x, y)
                    or self.is_oriented_by_meek_rule_2(x, y)
                    or self.is_oriented_by_meek_rule_3(x, y)
                ):
                    return False
            for y in self.neighbors[x]:
                if (
                    self.is_part_of_v_structure(x, y)
                    or self.is_oriented_by_meek_rule_1(x, y)
                    or self.is_oriented_by_meek_rule_2(x, y)
                    or self.is_oriented_by_meek_rule_3(x, y)
                ):
                    return False
                if (
                    self.is_part_of_v_structure(y, x)
                    or self.is_oriented_by_meek_rule_1(y, x)
                    or self.is_oriented_by_meek_rule_2(y, x)
                    or self.is_oriented_by_meek_rule_3(y, x)
                ):
                    return False
        return True

    def get_dag_extension(self) -> "PDAG":
        """
        Get a DAG in the Markov equivalence class of the PDAG.

        The returned DAG is also a PDAG object, but with only directed edges.

        Returns
        -------
        PDAG
            A DAG in the Markov equivalence class of the PDAG.
        """
        dag_extension = copy.deepcopy(self)
        dag_tmp = copy.deepcopy(self)
        nodes_tmp = set(self.nodes)

        while nodes_tmp:
            x = -1
            for node in nodes_tmp:
                if not dag_tmp.children[node]:
                    is_dag_extension = True
                    for neighbor in dag_tmp.neighbors[node]:
                        adjacent_neighbor = dag_tmp.adjacent[neighbor].copy()
                        adjacent_neighbor.add(neighbor)
                        if not dag_tmp.adjacent[node].issubset(adjacent_neighbor):
                            is_dag_extension = False
                            break
                    if is_dag_extension:
                        x = node
                        break
            if x == -1:
                raise RuntimeError("No consistent extension possible")
            while dag_tmp.neighbors[x]:
                neighbor = next(iter(dag_tmp.neighbors[x]))
                dag_tmp.remove_undirected_edge(neighbor, x)
                dag_extension.remove_undirected_edge(neighbor, x)
                dag_extension.add_directed_edge(neighbor, x)
            while dag_tmp.parents[x]:
                parent = next(iter(dag_tmp.parents[x]))
                dag_tmp.remove_directed_edge(parent, x)
            nodes_tmp.remove(x)
        return dag_extension

    def to_networkx(self):
        """
        Convert the PDAG to a networkx DiGraph object, where the PDAG undirected edges are
        represented using two directed edges in opposite directions.
        Call `get_dag_extension().to_networkx()` to get a networkx DiGraph object representing one
        of the DAG in the Markov equivalence class of the PDAG.

        Returns
        -------
        nx.DiGraph
            The networkx DiGraph object representing the PDAG.

        See Also
        --------
        get_dag_extension : Get a DAG in the Markov equivalence class of the PDAG.
        """
        import networkx as nx

        G = nx.DiGraph()
        for x in range(self.num_nodes):
            G.add_node(x)
        for x in range(self.num_nodes):
            for y in self.children[x]:
                G.add_edge(x, y)
        for x in range(self.num_nodes):
            for y in self.neighbors[x]:
                G.add_edge(x, y)
        return G

    def to_adjacency_matrix(self):
        """
        Convert the PDAG to an adjacency matrix.

        A directed edge (x, y) is represented by a 1 in cell (x, y) and a 0 in cell (y, x).
        An undirected edge {x, y} is represented by a 1 in cell (x, y) and a 1 in cell (y, x).
        Other cells are 0.

        Call `get_dag_extension().to_adjacency_matrix()` to get the adjacency matrix of a DAG
        in the Markov equivalence class of the PDAG.

        Returns
        -------
        numpy.ndarray
            The adjacency matrix of the PDAG.

        See Also
        --------
        get_dag_extension : Get a DAG in the Markov equivalence class of the PDAG.
        """
        adjacency_matrix = np.zeros((self.num_nodes, self.num_nodes))
        for x in range(self.num_nodes):
            for y in self.children[x]:
                adjacency_matrix[x, y] = 1
        for x in range(self.num_nodes):
            for y in self.neighbors[x]:
                adjacency_matrix[x, y] = 1
                adjacency_matrix[y, x] = 1
        return adjacency_matrix

    def __str__(self):
        undirected_edges_str = ", ".join(
            f"{x} - {y}" for x in range(self.num_nodes) for y in self.neighbors[x] if x < y
        )
        directed_edges_str = ", ".join(
            f"{x} → {y}" for x in range(self.num_nodes) for y in self.children[x]
        )
        return (
            f"PDAG: "
            f"undirected edges = {{{undirected_edges_str}}}, "
            f"directed edges = {{{directed_edges_str}}}"
        )

    def __repr__(self):
        return str(self)

    @staticmethod
    def from_networkx(G):
        """
        Create a PDAG from a networkx DiGraph object.

        The networkx DiGraph object should represent the PDAG, where the PDAG undirected edges are
        represented using two directed edges in opposite directions.

        Parameters
        ----------
        G : nx.DiGraph
            The networkx DiGraph object representing the PDAG.

        Returns
        -------
        PDAG
            The PDAG object representing the networkx DiGraph.

        See Also
        --------
        to_networkx : Convert the PDAG to a networkx DiGraph object.
        """
        pdag = PDAG(G.number_of_nodes())
        for x, y in G.edges:
            pdag.add_directed_edge(x, y)
        pdag.complete_cpdag()
        return pdag

    def complete_cpdag(self):
        pass
