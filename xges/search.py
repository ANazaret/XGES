import copy
import logging
import time
from collections import defaultdict
from typing import Dict, Set, Tuple

import numpy as np
from sortedcontainers import SortedListWithKey

from xges.bic_scorer import BICScorer
from xges.edge_queue_set import EdgeModificationsMap
from xges.operators import Delete, Insert, Reverse
from xges.pdag import PDAG
from xges.scorer import ScorerInterface

UnblockedPathsMap = Dict[Tuple[int, int], Set[Tuple[int, int]]]


class XGES:
    """
    The Extremely Greedy Equivalence Search algorithm for structure learning.

    Parameters
    ----------
    alpha: float (default=2.0)
        The penalty term for the BIC score.

    Attributes
    ----------
    alpha: float
        The penalty term for the BIC score.
    n_variables: int
        The number of variables in the data.
    scorer: Scorer
        The scorer used to evaluate the local scores. Default is BICScorerFast.
    pdag: PDAG
        The PDAG object representing the Markov equivalence class of the graph.
    initial_score: float
        The score of the empty graph.
    total_score: float
        The score of the current (final) graph.
    """

    def __init__(self, alpha: float = 2.0):
        self.alpha = alpha
        self.n_variables = None
        self.n_observations = None
        self.scorer = None
        self.pdag = None
        self.initial_score = None
        self.total_score = None
        self.statistics = None

        logger = logging.getLogger("xges")
        # set a default handler to console if not already set
        if not logger.handlers:
            logger.addHandler(logging.StreamHandler())
            logger.handlers[0].setFormatter(logging.Formatter("XGES: %(message)s"))
        self._logger = logger

    def fit(
        self,
        X: np.ndarray,
        extended_search: bool = True,
        use_fast_numba: bool = True,
        scorer: ScorerInterface = None,
        verbose: int = 1,
    ) -> PDAG:
        """
        Fit the XGES algorithm to the data.

        Parameters
        ----------
        X: np.ndarray
            The data matrix, where rows are samples and columns are variables.
        extended_search: bool (default=True)
            Whether to perform the extended search (XGES) or not (XGES-0).
        use_fast_numba: bool
            Whether to use the fast Numba implementation of the BIC scorer.
            Ignored if scorer is provided.
        scorer: ScorerInterface
            The scorer object to use. If None, the BIC scorer is used.
        verbose: int
            The verbosity level.
            0 is silent, 1 provides information before and after the search,
            2 provides information during the search and 3 provides very detailed information
            during the search. Default is 1.

        Returns
        -------
        PDAG:
            The PDAG object representing the Markov equivalence class of the graph.
        """
        if verbose == 0:
            self._logger.setLevel(logging.CRITICAL)
        elif verbose == 1:
            self._logger.setLevel(logging.INFO)
        elif verbose == 2:
            self._logger.setLevel(logging.DEBUG)
        elif verbose == 3:
            self._logger.setLevel(5)

        self._initialize_from_data(X, use_fast_numba=use_fast_numba, scorer=scorer)
        return self.fit_xges(extended_search)
        res = self.fit_xges(extended_search)
        self.graph_fitted = True
        return res

    def get_pdag(self) -> PDAG:
        """
        Get the PDAG object representing the Markov equivalence class of the graph.

        Returns
        -------
        PDAG
            The PDAG object.
        """
        return copy.deepcopy(self.pdag)

    def get_a_dag(self):
        """
        Get one of the DAGs in the Markov equivalence class of the graph.
        """
        dag_extension = self.pdag.get_dag_extension()
        return dag_extension

    def get_a_dag_networkx(self):
        """
        Get one of the DAGs in the Markov equivalence class of the graph, as a networkx object.

        """
        dag_extension = self.pdag.get_dag_extension()
        return dag_extension.to_networkx()

    def get_adjacency_matrix(self):
        """
        Get the adjacency matrix of the graph.
        """
        return self.pdag.to_adjacency_matrix()

    def _initialize_from_data(
        self,
        data: np.array,
        use_fast_numba: bool = True,
        scorer: ScorerInterface = None,
    ):
        """
        Initialize the XGES class from the data.

        See fit() for the parameters.
        """
        self.n_observations = data.shape[0]
        self.n_variables = data.shape[1]
        self.pdag = PDAG(self.n_variables)
        if scorer is not None:
            self.scorer = scorer
        elif use_fast_numba:
            try:
                from xges.bic_scorer_fast import BICScorerFast

                self.scorer = BICScorerFast(data, alpha=self.alpha)
            except ImportError:
                self._logger.warning(
                    "Package `numba` not found. "
                    "Falling back to the slower BICScorer implementation. "
                    "Install `numba` for better performance with `pip install numba`."
                )
                self.scorer = BICScorer(data, alpha=self.alpha)
        else:
            self.scorer = BICScorer(data, alpha=self.alpha)
        self.initial_score = self.scorer.score_pdag(self.pdag)
        self.total_score = self.initial_score
        self.statistics = defaultdict(int)

    def fit_xges(self, extended_search: bool = True) -> PDAG:
        """
        Fit the XGES (or XGES-0) algorithm to the data.

        Parameters
        ----------
        extended_search: bool (default=True)
            Whether to perform the extended search (XGES) or not (XGES-0).
        """
        candidate_inserts = SortedListWithKey(key=lambda op: -op.score)
        candidate_reverses = SortedListWithKey(key=lambda op: -op.score)
        candidate_deletes = SortedListWithKey(key=lambda op: -op.score)
        unblocked_paths_map = defaultdict(set)

        self._heuristic_xges0(
            candidate_inserts,
            candidate_reverses,
            candidate_deletes,
            unblocked_paths_map,
            initialize_inserts=True,
        )
        if extended_search:
            self._block_each_edge_and_search(unblocked_paths_map)

        self._logger.info(f"Final score: {self.total_score}")

        return self.get_pdag()

    def _heuristic_xges0(
        self,
        candidate_inserts: SortedListWithKey,
        candidate_reverses: SortedListWithKey,
        candidate_deletes: SortedListWithKey,
        unblocked_paths_map: UnblockedPathsMap,
        initialize_inserts: bool,
    ):
        """
        The XGES-0 algorithm.

        Parameters
        ----------
        candidate_inserts: SortedListWithKey
            The list of candidate inserts (usually empty at the beginning).
        candidate_reverses: SortedListWithKey
            The list of candidate reverses (usually empty at the beginning).
        candidate_deletes: SortedListWithKey
            The list of candidate deletes (usually empty at the beginning).
        unblocked_paths_map: UnblockedPathsMap
            The map of unblocked paths where the key is a tuple (x, y) and the value is a
            set of tuples (x, y) representing one unblocked path found from x to y.
        initialize_inserts: bool
            Whether to initialize the inserts or not.
            If True, the candidate_inserts list is filled with the possible inserts.
        """
        if initialize_inserts:
            if self.pdag.is_empty():
                self._find_inserts_on_empty_graph(candidate_inserts)
            else:
                # find all possible inserts
                for y in range(self.n_variables):
                    self._find_inserts_to_y(y, candidate_inserts, None, True)
        edge_modifications = EdgeModificationsMap()
        i_operations = 1

        last_insert = Insert(-1, -1, set(), -1, set())

        # XGES-0 main loop, in order: delete, reverse, insert; one operator per iteration
        while candidate_inserts or candidate_reverses or candidate_deletes:
            edge_modifications.clear()

            if candidate_deletes:
                # apply the best delete if possible
                best_delete = candidate_deletes.pop(0)
                if self.pdag.is_delete_valid(best_delete):
                    self.pdag.apply_delete(best_delete, edge_modifications)
                    self.total_score += best_delete.score
                    self._logger.debug(f"{i_operations}: {best_delete}")
                else:
                    continue
            elif candidate_reverses:
                # apply the best reverse if possible (no delete available)
                best_reverse = candidate_reverses.pop(0)
                if self.pdag.is_reverse_valid(best_reverse, unblocked_paths_map):
                    self.pdag.apply_reverse(best_reverse, edge_modifications)
                    self.total_score += best_reverse.score
                    self._logger.debug(f"{i_operations}: {best_reverse}")
                else:
                    continue
            elif candidate_inserts:
                # apply the best insert if possible (no delete or reverse available)
                best_insert = candidate_inserts.pop(0)
                if (
                    best_insert.y == last_insert.y
                    and abs(best_insert.score - last_insert.score) < 1e-10
                    and best_insert.x == last_insert.x
                    and best_insert.T == last_insert.T
                ):
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
            for edge_modification in edge_modifications:
                # log with level 5
                self._logger.log(5, f"Edge modification: {edge_modification}")

            # update the new possible operators
            self._update_operator_candidates_efficient(
                edge_modifications,
                candidate_inserts,
                candidate_reverses,
                candidate_deletes,
                unblocked_paths_map,
            )

    def _block_each_edge_and_search(self, unblocked_paths_map: UnblockedPathsMap):
        """
        Perform the extended search of the XGES algorithm.

        It applies each deletion and re-apply XGES-0 on the graph with the deletion enforced.

        Parameters
        ----------
        unblocked_paths_map: UnblockedPathsMap
            The map of unblocked paths obtained from the previous XGES-0 run.
        """
        all_edge_deletes = SortedListWithKey(key=lambda op: -op.score)
        deletes_of_pdag_are_updated = False
        index = 0

        while len(all_edge_deletes) or not deletes_of_pdag_are_updated:
            if index >= len(all_edge_deletes):
                all_edge_deletes.clear()
                for y in self.pdag.get_nodes():
                    self._find_deletes_to_y(y, all_edge_deletes, False)
                if len(all_edge_deletes) == 0:
                    break
                deletes_of_pdag_are_updated = True

            delete_ = all_edge_deletes.pop(0)
            if not self.pdag.is_delete_valid(delete_):
                continue
            # Apply the Delete
            xges_copy = copy.deepcopy(self)
            edge_modifications = EdgeModificationsMap()
            candidate_inserts = SortedListWithKey(key=lambda op: -op.score)
            candidate_reverses = SortedListWithKey(key=lambda op: -op.score)
            candidate_deletes = SortedListWithKey(key=lambda op: -op.score)
            xges_copy.pdag.apply_delete(delete_, edge_modifications)
            xges_copy.total_score += delete_.score
            xges_copy.pdag.add_forbidden_insert(delete_.x, delete_.y)
            blocked_paths_map_copy = copy.deepcopy(unblocked_paths_map)
            self._logger.debug(
                f"EXTENDED SEARCH {len(all_edge_deletes)} left (before updating deletes)"
            )
            for snd in edge_modifications:
                self._logger.log(5, f"\tEdge {snd}")

            xges_copy._update_operator_candidates_efficient(
                edge_modifications,
                candidate_inserts,
                candidate_reverses,
                candidate_deletes,
                blocked_paths_map_copy,
            )
            xges_copy._heuristic_xges0(
                candidate_inserts,
                candidate_reverses,
                candidate_deletes,
                blocked_paths_map_copy,
                initialize_inserts=False,
            )
            if self.pdag == xges_copy.pdag:
                continue
            if xges_copy.total_score - self.total_score > 1e-7 * self.n_observations:
                increase = xges_copy.total_score - self.total_score
                self._logger.debug(
                    f"EXTENDED SEARCH ACCEPTED: with increase {increase} and {delete_}"
                )
                self.total_score = xges_copy.total_score
                self.pdag = xges_copy.pdag
                unblocked_paths_map = blocked_paths_map_copy
                deletes_of_pdag_are_updated = False
                self.statistics["extended_search-accepted"] += 1
            else:
                self._logger.debug(
                    f"EXTENDED SEARCH REJECTED: {delete_} {xges_copy.total_score}"
                )
                self.statistics["extended_search-rejected"] += 1

    def _update_operator_candidates_naive(
        self,
        candidate_inserts: SortedListWithKey,
        candidate_reverses: SortedListWithKey,
        candidate_deletes: SortedListWithKey,
    ):
        """
        Update the candidate operators of all variables. This is the naive implementation used
        by GES.

        XGES uses a more efficient implementation in `_update_operator_candidates_efficient`.

        Parameters
        ----------
        candidate_inserts: SortedListWithKey
            The list of candidate inserts, which will be updated.
        candidate_reverses: SortedListWithKey
            The list of candidate reverses, which will be updated.
        candidate_deletes: SortedListWithKey
            The list of candidate deletes, which will be updated.

         See Also
        --------
        _update_operator_candidates_efficient
        """
        candidate_inserts.clear()
        candidate_reverses.clear()
        candidate_deletes.clear()
        for y in range(self.n_variables):
            self._find_inserts_to_y(y, candidate_inserts)
            self._find_reverse_to_y(y, candidate_reverses)
            self._find_deletes_to_y(y, candidate_deletes)

    def _update_operator_candidates_efficient(
        self,
        edge_modifications: EdgeModificationsMap,
        candidate_inserts,
        candidate_reverses,
        candidate_deletes,
        unblocked_paths_map: Dict[Tuple[int, int], Set[Tuple[int, int]]],
    ):
        start_time = time.time()

        # First, undo all the edge modifications
        for edge_modification in edge_modifications:
            self.pdag.apply_edge_modification(edge_modification, True)

        full_insert_to_y = set()
        partial_insert_to_y = defaultdict(set)
        full_delete_to_y = set()
        full_delete_from_x = set()
        delete_x_y = set()
        full_reverse_to_y = set()
        full_reverse_from_x = set()
        reverse_x_y = set()

        # Re-apply the edge modifications one by one and update the operators
        for edge_modification in edge_modifications:

            if edge_modification.is_old_directed():
                a = edge_modification.get_old_source()
                b = edge_modification.get_old_target()
            elif edge_modification.is_new_directed():
                a = edge_modification.get_new_source()
                b = edge_modification.get_new_target()
            else:
                a = edge_modification.x
                b = edge_modification.y

            # Track inserts
            modification_id = edge_modification.get_modification_id()
            if modification_id in [1, 2]:
                # 1. a  b becomes a -- b
                # 2. a  b becomes a → b
                if modification_id == 1:
                    # y = a
                    full_insert_to_y.add(a)
                # y = b
                full_insert_to_y.add(b)
                # y \in Ne(a) ∩ Ne(b)
                full_insert_to_y.update(
                    self.pdag.get_neighbors(a) & self.pdag.get_neighbors(b)
                )
                # x=a and y \in Ne(b)
                for target in self.pdag.get_neighbors(b):
                    partial_insert_to_y[target].add(a)
                # x=b and y \in Ne(a)
                for target in self.pdag.get_neighbors(a):
                    partial_insert_to_y[target].add(b)
            elif modification_id == 3:
                # a -- b becomes a  b
                # x=a and y \in Ne(b) u {b}
                for target in self.pdag.get_neighbors(b):
                    if target != a:
                        partial_insert_to_y[target].add(a)
                partial_insert_to_y[b].add(a)
                # y=a and x \in Ad(b)
                partial_insert_to_y[a].update(self.pdag.get_adjacent(b))
                # x=b and y \in Ne(a) u {a}
                for target in self.pdag.get_neighbors(a):
                    if target != b:
                        partial_insert_to_y[target].add(b)
                partial_insert_to_y[a].add(b)
                # y=b and x \in Ad(a)
                partial_insert_to_y[b].update(self.pdag.get_adjacent(a))
                # SD(x,y,a,b)
                if (a, b) in unblocked_paths_map:
                    for x, y in unblocked_paths_map[(a, b)]:
                        partial_insert_to_y[y].add(x)
                # SD(x,y,b,a)
                if (b, a) in unblocked_paths_map:
                    for x, y in unblocked_paths_map[(b, a)]:
                        partial_insert_to_y[y].add(x)
            elif modification_id == 4:
                # a -- b becomes a → b
                # y = a and x \in Ad(b)
                partial_insert_to_y[a].update(self.pdag.get_adjacent(b))
                # y = b
                full_insert_to_y.add(b)
                # SD(x,y,b,a)
                if (b, a) in unblocked_paths_map:
                    for x, y in unblocked_paths_map[(b, a)]:
                        partial_insert_to_y[y].add(x)
            elif modification_id == 5:
                # a → b becomes a  b
                # x=a and y \in Ne(b) u {b}
                for target in self.pdag.get_neighbors(b):
                    partial_insert_to_y[target].add(a)
                partial_insert_to_y[b].add(a)
                # x=b and y \in Ne(a) u {a}
                for target in self.pdag.get_neighbors(a):
                    partial_insert_to_y[target].add(b)
                partial_insert_to_y[a].add(b)
                full_insert_to_y.add(b)
                # SD(x,y,a,b)
                if (a, b) in unblocked_paths_map:
                    for x, y in unblocked_paths_map[(a, b)]:
                        partial_insert_to_y[y].add(x)
            elif modification_id == 6:
                # a → b becomes a -- b
                full_insert_to_y.add(a)
                full_insert_to_y.add(b)
            elif modification_id == 7:
                # a → b becomes a ← b
                full_insert_to_y.add(a)
                full_insert_to_y.add(b)
                # SD(x,y,a,b)
                if (a, b) in unblocked_paths_map:
                    for x, y in unblocked_paths_map[(a, b)]:
                        partial_insert_to_y[y].add(x)

            # Track deletes
            if modification_id in [1, 2]:
                # a  b becomes a -- b
                if modification_id == 1:
                    full_delete_to_y.add(a)
                full_delete_to_y.add(b)
                full_delete_from_x.add(a)
                full_delete_from_x.add(b)

                # x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b) [almost never happens]
                x_intersection = self.pdag.get_adjacent(a) & self.pdag.get_adjacent(b)
                if x_intersection:
                    y_intersection = self.pdag.get_neighbors(
                        a
                    ) & self.pdag.get_neighbors(b)
                    delete_x_y.update(
                        (x, y) for x in x_intersection for y in y_intersection
                    )
            elif modification_id in [4, 5]:
                # a -- b becomes a → b
                # a → b becomes a  b
                full_delete_to_y.add(b)
            elif modification_id in [6, 7]:
                # a → b becomes a -- b
                # a → b becomes a ← b
                full_delete_to_y.add(a)
                full_delete_to_y.add(b)

            # Track reverse
            if modification_id in [1, 2]:
                # 1. a  b becomes a -- b
                # 2. a  b becomes a → b
                # y \in {a, b}
                if modification_id == 1:
                    full_reverse_to_y.add(a)
                full_reverse_to_y.add(b)
                # y \in Ne(a) ∩ Ne(b)
                full_reverse_to_y.update(
                    self.pdag.get_neighbors(a) & self.pdag.get_neighbors(b)
                )
                # x \in {a, b}
                full_reverse_from_x.add(a)
                full_reverse_from_x.add(b)
            elif modification_id == 3:
                # a -- b becomes a  b
                # y \in {a, b} or x \in {a, b}
                full_reverse_to_y.add(a)
                full_reverse_to_y.add(b)
                full_reverse_from_x.add(a)
                full_reverse_from_x.add(b)
                # SD(x,y,a,b) or SD(x,y,b,a)
                for path_key in [(a, b), (b, a)]:
                    if path_key in unblocked_paths_map:
                        reverse_x_y.update(unblocked_paths_map[path_key])
                        unblocked_paths_map.pop(path_key)
            elif modification_id == 4:
                # a -- b becomes a → b
                # y \in {a, b} or x = b
                full_reverse_to_y.add(a)
                full_reverse_to_y.add(b)
                full_reverse_from_x.add(b)
                # SD(x,y,b,a)
                if (b, a) in unblocked_paths_map:
                    reverse_x_y.update(unblocked_paths_map[(b, a)])
                    unblocked_paths_map.pop((b, a))
            elif modification_id == 5:
                # a → b becomes a  b
                # y = b or x \in {a, b}
                full_reverse_to_y.add(b)
                full_reverse_from_x.add(a)
                full_reverse_from_x.add(b)
                # SD(x,y,a,b)
                if (a, b) in unblocked_paths_map:
                    reverse_x_y.update(unblocked_paths_map[(a, b)])
                    unblocked_paths_map.pop((a, b))
            elif modification_id == 6:
                # a → b becomes a -- b
                # y \in {a, b} or x = b
                full_reverse_to_y.add(a)
                full_reverse_to_y.add(b)
                full_reverse_from_x.add(b)
            elif modification_id == 7:
                # a → b becomes a ← b
                # y \in {a, b} or x \in {a, b}
                full_reverse_to_y.add(a)
                full_reverse_to_y.add(b)
                full_reverse_from_x.add(a)
                full_reverse_from_x.add(b)
                # SD(x,y,a,b)
                if (a, b) in unblocked_paths_map:
                    reverse_x_y.update(unblocked_paths_map[(a, b)])
                    unblocked_paths_map.pop((a, b))

            self.pdag.apply_edge_modification(edge_modification)

        # Find the inserts
        # Step 1: remove partial inserts that are now full
        keys_to_erase = [y for y in partial_insert_to_y if y in full_insert_to_y]
        for key in keys_to_erase:
            partial_insert_to_y.pop(key)

        # Step 2: find the partial inserts
        for y, xs in partial_insert_to_y.items():
            for x in xs:
                if x not in self.pdag.get_adjacent(y) and x != y:
                    # self._logger.debug(f"Partial insert: x={x}, y={y}")
                    self._find_inserts_to_y(y, candidate_inserts, x)

        # Step 3: find the full inserts
        for y in full_insert_to_y:
            # self._logger.debug(f"Full insert: y={y}")
            self._find_inserts_to_y(y, candidate_inserts)

        # Find the deletes
        # Step 1: form the edges to delete
        for x in full_delete_from_x:
            delete_x_y.update((x, y) for y in self.pdag.get_neighbors(x))
            delete_x_y.update((x, y) for y in self.pdag.get_children(x))
        for y in full_delete_to_y:
            delete_x_y.update((x, y) for x in self.pdag.get_parents(y))
            delete_x_y.update((x, y) for x in self.pdag.get_neighbors(y))

        # Step 2: find the deletes
        for x, y in delete_x_y:
            if x != y:
                self._find_delete_to_y_from_x(y, x, candidate_deletes)

        # Find the reverses
        # Step 1: form the edges to reverse
        for x in full_reverse_from_x:
            reverse_x_y.update((x, y) for y in self.pdag.get_parents(x))
        for y in full_reverse_to_y:
            reverse_x_y.update((x, y) for x in self.pdag.get_children(y))

        # Step 2: find the reverses
        for x, y in reverse_x_y:
            # check that x ← y
            if self.pdag.has_directed_edge(y, x) and x != y:
                self._find_reverse_to_y_from_x(y, x, candidate_reverses)

        self.statistics["update_operators-time"] += time.time() - start_time

    def _find_inserts_to_y(
        self, y, candidate_inserts, parent_x=None, positive_only=True
    ):
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
            neighbors_y_not_adjacent_x = list(
                self.pdag.get_neighbors_not_adjacent(y, x)
            )

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
                    if T.issubset(adjacent_z) and neighbors_y_adjacent_x.issubset(
                        adjacent_z
                    ):
                        T_prime = T.copy()
                        T_prime.add(z)
                        effective_parents_prime = effective_parents.copy()
                        effective_parents_prime.add(z)
                        stack.append((T_prime, idx + i + 1, effective_parents_prime))

    def _find_delete_to_y_from_x(self, y, x, candidate_deletes, positive_only=True):
        parents_y = self.pdag.get_parents(y)
        neighbors_y_adjacent_x = list(self.pdag.get_neighbors_adjacent(y, x))
        directed_xy = self.pdag.has_directed_edge(x, y)

        stack = [(set(), 0, parents_y.union({x}))]

        while stack:
            C, idx, effective_parents = stack.pop()

            score = self.scorer.score_delete(y, list(effective_parents), x)
            if score > 0 or not positive_only:
                candidate_deletes.add(
                    Delete(x, y, C, score, effective_parents, directed_xy)
                )

            for i, z in enumerate(neighbors_y_adjacent_x[idx:]):
                adjacent_z = self.pdag.get_adjacent(z)
                if C.issubset(adjacent_z):
                    C_prime = C.copy()
                    C_prime.add(z)
                    effective_parents_prime = effective_parents.copy()
                    effective_parents_prime.add(z)
                    stack.append((C_prime, idx + 1 + i, effective_parents_prime))

    def _find_deletes_to_y(self, y, candidate_deletes, positive_only=True):
        neighbors_y = self.pdag.get_neighbors(y)
        parents_y = self.pdag.get_parents(y)

        for x in parents_y:
            self._find_delete_to_y_from_x(y, x, candidate_deletes, positive_only)

        for x in neighbors_y:
            self._find_delete_to_y_from_x(y, x, candidate_deletes, positive_only)

    def _find_reverse_to_y(self, y, candidate_reverses):
        # look for all possible x ← y
        children_y = self.pdag.get_children(y)

        for x in children_y:
            parents_x = self.pdag.get_parents(x)
            candidate_inserts = set()
            self._find_inserts_to_y(y, candidate_inserts, x, positive_only=False)

            for insert in candidate_inserts:
                score = insert.score + self.scorer.score_delete(x, parents_x, y)
                if score > 0:
                    candidate_reverses.add(Reverse(insert, score, parents_x.copy()))

    def _find_reverse_to_y_from_x(self, y, x, candidate_reverses):
        if not self.pdag.has_directed_edge(y, x):
            return
        candidate_inserts = set()
        self._find_inserts_to_y(y, candidate_inserts, x, positive_only=False)
        parents_x = self.pdag.get_parents(x)
        for insert in candidate_inserts:
            score = insert.score + self.scorer.score_delete(x, parents_x, y)
            if score > 0:
                candidate_reverses.add(Reverse(insert, score, parents_x.copy()))

    def _find_inserts_on_empty_graph(self, candidate_inserts):
        paired_scores = self.scorer.score_all_pairs()
        for x in range(self.n_variables):
            for y in range(self.n_variables):
                if x == y:
                    continue
                # we consider both x → y and y → x
                score = paired_scores[x, y]
                if score > 0:
                    candidate_inserts.add(Insert(x, y, set(), score, set()))
