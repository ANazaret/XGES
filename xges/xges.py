"""
Adaptation from c++ of:

class XGES {
public:
    XGES(int n_variables, ScorerInterface *scorer);
    XGES(const XGES &other);

    void fit_xges(bool extended_search);
    void fit_ops(bool use_reverse);
    void fit_ges(bool use_reverse);
    double get_score() const;
    double get_initial_score() const;
    const PDAG &get_pdag() const;

    std::unique_ptr<PDAG> ground_truth_pdag;
    std::map<std::string, double> statistics;

private:

    void heuristic_xges0(std::vector<Insert> &candidate_inserts,
                         std::vector<Reverse> &candidate_reverses,
                         std::vector<Delete> &candidate_deletes,
                         UnblockedPathsMap &unblocked_paths_map,
                         bool initialize_inserts = true);
    void update_operator_candidates_naive(std::vector<Insert> &candidate_inserts,
                                          std::vector<Reverse> &candidate_reverses,
                                          std::vector<Delete> &candidate_deletes);
    void update_operator_candidates_efficient(EdgeModificationsMap &edge_modifications,
                                              std::vector<Insert> &candidate_inserts,
                                              std::vector<Reverse> &candidate_reverses,
                                              std::vector<Delete> &candidate_deletes,
                                              UnblockedPathsMap &unblocked_paths_map);
    void block_each_edge_and_research(UnblockedPathsMap &unblocked_paths_map);

};

"""
from xges.pdag import PDAG
from xges.operators import Insert, Delete, Reverse


class XGES:
    def __init__(self, n_variables, scorer):
        self.n_variables = n_variables
        self.scorer = scorer
        self.pdag = PDAG(n_variables)
        self.initial_score = ...
        self.total_score = self.initial_score
        self._logger = None

    def fit(self, X):
        return self.fit_xges()

    def fit_xges(self, X, extended_search=True):
        """
        Fit the XGES algorithm to the data.
        """
        self._heuristic_xges0(X, initialize_inserts=True)
        if extended_search:
            self._block_each_edge_and_search()

    def _heuristic_xges0(self, candidate_inserts, candidate_reverses, candidate_deletes, unblocked_paths_map,
                         initialize_inserts):
        """

        """
        if initialize_inserts:
            # TODO plug the initialization on empty PDAG
            # find all possible inserts
            for y in range(self.n_variables):
                self.find_inserts_to_y(y, candidate_inserts, None, True)

        edge_modifications = {}
        i_operations = 1

        last_insert = Insert(-1, -1, set(), -1, set())

        # XGES-0 main loop, in order: delete, reverse, insert; one operator per iteration
        while candidate_inserts or candidate_reverses or candidate_deletes:
            edge_modifications.clear()

            if candidate_deletes:
                # apply the best delete if possible
                best_delete = candidate_deletes.pop()
                if self.pdag.is_delete_valid(best_delete):
                    self.pdag.apply_delete(best_delete, edge_modifications)
                    self.total_score += best_delete.score
                    self._logger.debug(f"{i_operations}: {best_delete}")
                else:
                    continue
            elif candidate_reverses:
                # apply the best reverse if possible (no delete available)
                best_reverse = candidate_reverses.pop()
                if self.pdag.is_reverse_valid(best_reverse, unblocked_paths_map):
                    self.pdag.apply_reverse(best_reverse, edge_modifications)
                    self.total_score += best_reverse.score
                    self._logger.debug(f"{i_operations}: {best_reverse}")
                else:
                    continue
            elif candidate_inserts:
                # apply the best insert if possible (no delete or reverse available)
                best_insert = candidate_inserts.pop()
                if best_insert.y == last_insert.y and abs(
                        best_insert.score - last_insert.score) < 1e-10 and best_insert.x == last_insert.x and best_insert.T == last_insert.T:
                    self.statistics["probable_insert_duplicates"] += 1
                    continue
                last_insert = best_insert
                if self.pdag.is_insert_valid(last_insert, unblocked_paths_map):
                    self.pdag.apply_insert(last_insert, edge_modifications)
                    self.total_score += last_insert.score
                    self._logger.debug(f"{i_operations}: {last_insert}")
                else:
                    continue

            # if we reach this point, we have applied an operator
            i_operations += 1
            for edge_modification in edge_modifications.values():
                self._logger.trace(f"\tEdge {edge_modification}")

            # update the new possible operators
            # self.update_operator_candidates_efficient(edge_modifications, candidate_inserts, candidate_reverses,
            #                                           candidate_deletes, unblocked_paths_map)
            self.update_operator_candidates_naive(candidate_inserts, candidate_reverses, candidate_deletes)

    def update_operator_candidates_naive(self, candidate_inserts, candidate_reverses, candidate_deletes):
        candidate_inserts.clear()
        candidate_reverses.clear()
        candidate_deletes.clear()
        for y in range(self.n_variables):
            self.find_inserts_to_y(y, candidate_inserts)
            self.find_reverse_to_y(y, candidate_reverses)
            self.find_deletes_to_y(y, candidate_deletes)

    def find_inserts_to_y(self, y, candidate_inserts, parent_x=None, positive_only=True):
        adjacent_y = self.pdag.get_adjacent(y)
        parents_y = self.pdag.get_parents(y)

        possible_parents = set()

        if parent_x is not None:
            possible_parents.add(parent_x)
        else:
            nodes = self.pdag.get_nodes()
            # 1. x is not adjacent to y (x ∉ Ad(y))
            possible_parents = nodes - adjacent_y
            possible_parents.remove(y)

        for x in possible_parents:
            neighbors_y_adjacent_x = self.pdag.get_neighbors_adjacent(y, x)
            # 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
            # So in particular, [Ne(y) ∩ Ad(x)] is a clique
            if not self.pdag.is_clique(neighbors_y_adjacent_x):
                continue

            # 2. T ⊆ Ne(y) \ Ad(x)
            neighbors_y_not_adjacent_x = list(self.pdag.get_neighbors_not_adjacent(y, x))

            effective_parents_y = neighbors_y_adjacent_x
            effective_parents_y.update(parents_y)

            stack = [(set(), 0, effective_parents_y)]

            while stack:
                T, idx, effective_parents = stack.pop()

                score = self.scorer.score_insert(y, effective_parents, x)
                if score > 0 or not positive_only:
                    candidate_inserts.add(Insert(x, y, T, score, effective_parents))

                for i, z in enumerate(neighbors_y_not_adjacent_x[idx:]):
                    adjacent_z = self.pdag.get_adjacent(z)
                    if T.issubset(adjacent_z) and neighbors_y_adjacent_x.issubset(adjacent_z):
                        T_prime = T.copy()
                        T_prime.add(z)
                        effective_parents_prime = effective_parents.copy()
                        effective_parents_prime.add(z)
                        stack.append((T_prime, idx + i + 1, effective_parents_prime))

    def find_delete_to_y_from_x(self, y, x, candidate_deletes, positive_only=True):
        parents_y = self.pdag.get_parents(y)
        neighbors_y_adjacent_x = self.pdag.get_neighbors_adjacent(y, x)
        directed_xy = self.pdag.has_directed_edge(x, y)

        stack = [(set(), 0, parents_y.union({x}))]

        while stack:
            C, idx, effective_parents = stack.pop()

            score = self.scorer.score_delete(y, list(effective_parents), x)
            if score > 0 or not positive_only:
                candidate_deletes.add(Delete(x, y, C, score, effective_parents, directed_xy))

            for i, z in enumerate(neighbors_y_adjacent_x[idx:]):
                adjacent_z = self.pdag.get_adjacent(z)
                if C.issubset(adjacent_z):
                    C_prime = C.copy()
                    C_prime.add(z)
                    effective_parents_prime = effective_parents.copy()
                    effective_parents_prime.add(z)
                    stack.append((C_prime, idx + 1 + i, effective_parents_prime))

    def find_deletes_to_y(self, y, candidate_deletes, positive_only=True):
        neighbors_y = self.pdag.get_neighbors(y)
        parents_y = self.pdag.get_parents(y)

        for x in parents_y:
            self.find_delete_to_y_from_x(y, x, candidate_deletes, positive_only)

        for x in neighbors_y:
            self.find_delete_to_y_from_x(y, x, candidate_deletes, positive_only)

    def find_reverse_to_y(self, y, candidate_reverses):
        # look for all possible x ← y
        children_y = self.pdag.get_children(y)

        for x in children_y:
            parents_x = self.pdag.get_parents(x)
            candidate_inserts = []
            self.find_inserts_to_y(y, candidate_inserts, x, False)

            for insert in candidate_inserts:
                score = insert.score + self.scorer.score_delete(x, parents_x, y)
                if score > 0:
                    candidate_reverses.add(Reverse(insert, score, parents_x))

    def find_reverse_to_y_from_x(self, y, x, candidate_reverses):
        if not self.pdag.has_directed_edge(y, x):
            return
        candidate_inserts = []
        self.find_inserts_to_y(y, candidate_inserts, x, False)
        parents_x = self.pdag.get_parents(x)
        for insert in candidate_inserts:
            score = insert.score + self.scorer.score_delete(x, parents_x, y)
            if score > 0:
                candidate_reverses.add(Reverse(insert, score, parents_x))
