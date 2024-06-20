from collections import defaultdict

import numba as nb
import numpy as np
from scipy.special import gammaln

from xges.scorer import ScorerInterface


@nb.jit(
    nb.float64[:, ::1](
        nb.float64[:, ::1],  # the 1 means it is contiguous
        nb.int64,
        nb.float64,
    ),
    nopython=True,
    fastmath=True,
    parallel=True,
)
def score_all_pairs(cov, n, alpha=1.0):
    """
    Return a matrix of all the diff scores for all pairs of variables.
    In [i,j] is the score of inserting i in the parents of j.
    Insert(i,j) = score(j, [i]) - score(j, [])
                =`-n/2 log(1-Cov_{i,j}^2/(Cov_{i,i}Cov_{j,j}) - alpha log n/2`
    """
    diag_sqrt = np.sqrt(np.diag(cov))
    corr = cov / np.outer(diag_sqrt, diag_sqrt)
    scores = -0.5 * n * np.log(1 - corr**2) - 0.5 * np.log(n) * alpha
    return scores


@nb.jit(
    nb.float64[:, :](
        nb.float64[:, ::1],  # the 1 means it is contiguous
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


class BICScorerFast(ScorerInterface):
    def __init__(self, data, alpha):
        super().__init__()
        self.data = data
        self.alpha = alpha
        self.covariance_matrix = self.compute_covariance(data)
        self.n_variables = data.shape[1]
        self.n_samples = data.shape[0]
        self.cache = [defaultdict(float) for _ in range(self.n_variables)]
        self.statistics = defaultdict(int)

        self.cache_score_all_pairs = None

    @staticmethod
    def compute_covariance(data):
        n_samples = data.shape[0]
        centered = data - np.mean(data, axis=0)
        covariance_matrix = np.dot(centered.T, centered) / n_samples
        return covariance_matrix

    @staticmethod
    def log_binomial(n, k):
        log_n_choose_k = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        return log_n_choose_k

    def score_all_pairs(self):
        self.cache_score_all_pairs = score_all_pairs(
            self.covariance_matrix, self.n_samples, self.alpha
        )
        # todo: consider updating cache (it would require saving the intermediate scores in
        #  score_all_pairs)
        return self.cache_score_all_pairs

    def local_score(self, target, parents):
        self.statistics["local_score-#calls-total"] += 1
        parents = sorted(parents)
        cache_key = tuple(parents)
        if cache_key in self.cache[target]:
            return self.cache[target][cache_key]
        self.statistics["local_score-#calls-nocache"] += 1

        bic = local_bic_fast(
            target,
            np.array(parents, dtype=int),
            self.n_samples,
            self.covariance_matrix,
            self.alpha,
        )

        self.cache[target][cache_key] = bic

        return bic
