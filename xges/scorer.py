from abc import ABC, abstractmethod
from collections import defaultdict
from time import time

import numpy as np


class ScorerInterface(ABC):
    def __init__(self):
        self.statistics = defaultdict(int)

    @abstractmethod
    def local_score(self, target, parents):
        pass

    def score_all_pairs(self):
        scores = np.zeros((self.n_variables, self.n_variables))
        for i in range(self.n_variables):
            for j in range(self.n_variables):
                if i == j:
                    continue
                # i to j is the same as j to i
                parents = {i}
                scores[i, j] = self.local_score(j, parents) - self.local_score(j, set())

        return scores

    def score_insert(self, target, parents, parent_to_add):
        self.statistics["score_insert-#calls"] += 1
        start_time = time()
        score_without_new_parent = self.local_score(target, parents)
        parents_with_new_parent = parents.copy()
        parents_with_new_parent.add(parent_to_add)
        score_with_new_parent = self.local_score(target, parents_with_new_parent)

        self.statistics["score_insert-time"] += time() - start_time
        return score_with_new_parent - score_without_new_parent

    def score_delete(self, target, parents, parent_to_remove):
        self.statistics["score_delete-#calls"] += 1
        start_time = time()
        score_with_old_parent = self.local_score(target, parents)
        parents_without_old_parent = parents.copy()
        parents_without_old_parent.remove(parent_to_remove)
        score_without_old_parent = self.local_score(target, parents_without_old_parent)

        self.statistics["score_delete-time"] += time() - start_time
        return score_without_old_parent - score_with_old_parent

    def score_pdag(self, pdag):
        score = 0
        dag = pdag.get_dag_extension()
        # compute the score at nodes that are variables
        for target in dag.get_nodes():
            score += self.local_score(target, dag.get_parents(target))
        return score
