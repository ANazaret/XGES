import numpy as np
from scipy.special import gammaln
from collections import defaultdict
from time import time
from xges.scorer import ScorerInterface


class BICScorer(ScorerInterface):
    def __init__(self, data, alpha):
        super().__init__()
        self.data = data
        self.alpha = alpha
        self.covariance_matrix = self.compute_covariance(data)
        self.n_variables = data.shape[1]
        self.n_samples = data.shape[0]
        self.cache = [defaultdict(float) for _ in range(self.n_variables)]
        self.statistics = defaultdict(int)

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

    def local_score(self, target, parents):
        self.statistics["local_score-#calls-total"] += 1
        parents = sorted(parents)
        cache_key = tuple(parents)
        if cache_key in self.cache[target]:
            return self.cache[target][cache_key]
        self.statistics["local_score-#calls-nocache"] += 1

        cov_target_target = self.covariance_matrix[target, target]

        if len(parents) == 0:
            sigma = cov_target_target
        else:
            parents = np.array(parents)
            cov_parents_parents = self.covariance_matrix[np.ix_(parents, parents)]
            cov_parents_target = self.covariance_matrix[np.ix_(parents, [target])]

            beta = np.linalg.solve(cov_parents_parents, cov_parents_target)
            sigma = cov_target_target - np.dot(cov_parents_target.T, beta)

        log_likelihood_no_constant = -0.5 * self.n_samples * (1 + np.log(sigma))

        bic_regularization = 0.5 * np.log(self.n_samples) * (len(parents) + 1.) * self.alpha

        bic = log_likelihood_no_constant - bic_regularization

        self.cache[target][cache_key] = bic

        return bic
