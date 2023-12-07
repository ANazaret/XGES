import networkx as nx
import numba as nb
import numpy as np

from evaluation.pdag import PDAG


@nb.jit(
    nb.float64[:, :](
        nb.float64[:, ::1],
        nb.int64[::1],
        nb.int64[::1],
    ),
    nopython=True,
    fastmath=True,
)
def numba_ix(arr, rows, cols):
    """
    From: https://github.com/numba/numba/issues/5894#issuecomment-974701551

    Numba compatible implementation of arr[np.ix_(rows, cols)] for 2D arrays.

    Parameters
    ----------

    """
    one_d_index = np.zeros(len(rows) * len(cols), dtype=np.int32)
    for i, r in enumerate(rows):
        start = i * len(cols)
        one_d_index[start : start + len(cols)] = cols + arr.shape[1] * r

    arr_1d = arr.reshape((arr.shape[0] * arr.shape[1], 1))
    slice_1d = np.take(arr_1d, one_d_index)
    return slice_1d.reshape((len(rows), len(cols)))


nb.float64(nb.int64, nb.int64[:], nb.int64, nb.float64[:, :], nb.float64),


@nb.jit(
    nb.float64(
        nb.int64,
        nb.int64[::1],  # the 1 means it is contiguous
        nb.int64,
        nb.float64[:, ::1],  # the 1 means it is contiguous
        nb.float64,
    ),
    nopython=True,
    fastmath=True,
)
def local_bic_fast(target, parents, n, cov, alpha=1.0):
    cov_parents_parents = numba_ix(cov, parents, parents)
    cov_parents_target = cov[parents, target]
    cov_target_target = cov[target, target]
    beta = np.linalg.solve(cov_parents_parents, cov_parents_target)
    sigma = cov_target_target - cov_parents_target.T @ beta

    log_likelihood_no_constant = -0.5 * n * (1 + np.log(sigma))
    bic_regularization = 0.5 * np.log(n) * (len(parents) + 1) * alpha

    bic = log_likelihood_no_constant - bic_regularization
    return bic


class BICScorer:
    def __init__(
        self,
        n_samples,
        n_variables,
        alpha=1.0,
        data=None,
        intervention_labels=None,
        cov=None,
        ignored_variables=None,
    ):
        self.n, self.d = n_samples, n_variables
        self.alpha = alpha
        self.ignored_variables = ignored_variables or set()

        if cov is None:
            if intervention_labels is not None:
                # one-hot encode intervention labels (note that label -1 should give all zeros)
                intervention_labels = np.array(intervention_labels)
                intervention_labels = np.eye(intervention_labels.max() + 1)[intervention_labels]
                intervention_labels = intervention_labels[:, :-1]
                data = np.concatenate([data, intervention_labels], axis=1)
            cov = np.cov(data, rowvar=False, bias=True).reshape((data.shape[1], data.shape[1]))
            # reshape: if we have only one variable, cov is a scalar, but we want a 1x1 matrix
        self._cov = cov
        self._cache = dict()

    def check_cache(self, cache_key):
        # should sort parents to make sure cache is hit, but this is slow.
        # benchmark: if sorted, 51.0% hits vs unsorted is 50.6% hits ... so not worth it
        # cache_key = (target, tuple(sorted(parents)))
        if cache_key in self._cache:
            return self._cache[cache_key]
        return None

    def set_cache(self, cache_key, value):
        self._cache[cache_key] = value

    def local_diff_score(self, target, parents, candidate_parent):
        if target in self.ignored_variables:
            # return the bic score cost of adding 1 parent
            return -0.5 * np.log(self.n) * self.alpha

        parents_with_candidate = np.array(list(parents) + [candidate_parent])

        # score without candidate parent
        parents_only = parents_with_candidate[:-1]
        cache_key = (target, tuple(parents_only))
        base_score = self.check_cache(cache_key)
        if base_score is None:
            # base_score = self.local_score(target, parents_only)
            base_score = local_bic_fast(target, parents_only, self.n, self._cov, self.alpha)
            self.set_cache(cache_key, base_score)

        # score with candidate parent
        cache_key = (target, tuple(parents_with_candidate))
        candidate_score = self.check_cache(cache_key)
        if candidate_score is None:
            # candidate_score = self.local_score(target, parents_with_candidate)
            candidate_score = local_bic_fast(target, parents_with_candidate, self.n, self._cov, self.alpha)
            self.set_cache(cache_key, candidate_score)

        return candidate_score - base_score

    def score_of_dag(self, dag):
        """
        Return the score for a DAG.
        """
        if isinstance(dag, nx.DiGraph) and not nx.is_directed_acyclic_graph(dag):
            # it might be a valid PDAG
            dag = PDAG.from_digraph(dag)
        if isinstance(dag, PDAG):
            dag = dag.to_dag()

        score = 0
        for target in range(self.d):
            if target in self.ignored_variables:
                continue
            parents = list(dag.predecessors(target))
            local_score = local_bic_fast(target, np.array(parents, dtype=int), self.n, self._cov, self.alpha)
            score += local_score

        # score += self.d * (self.n + np.log(self.n)) / 2

        return score
