//
// Created by Achille Nazaret on 11/7/23.
//
#include "BICScorer.h"

MatrixXd compute_covariance(const MatrixXd &data) {
    const long n_samples = data.rows();

    const MatrixXd centered = data.rowwise() - data.colwise().mean();
    MatrixXd covariance_matrix = centered.adjoint() * centered / n_samples;
    return covariance_matrix;
}

MatrixXd compute_covariance(const MatrixXd &data, const VectorXi &interventions_index) {
    const long n_samples = data.rows();
    const MatrixXd centered = data.rowwise() - data.colwise().mean();
    MatrixXd covariance_matrix = centered.adjoint() * centered / n_samples;
    return covariance_matrix;
}

BICScorer::BICScorer(const MatrixXd &data, double alpha)
    : data(data), alpha(alpha), covariance_matrix(compute_covariance(data)) {
    n_variables = data.cols();
    n_samples = data.rows();
    cache.resize(n_variables);
}

double log_binomial(const int n, const int k) {
    // by definition: log(n choose k) = log(n!) - log(k!) - log((n-k)!)
    // and log-gamma(n+1) = log(n!)
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}


double BICScorer::local_score(int target, const FlatSet &parents) {
    // cache lookup
    statistics["local_score-#calls-total"]++;
    auto &cache_target = cache[target];
    auto it = cache_target.find(parents);
    if (it != cache_target.end()) { return it->second; }
    statistics["local_score-#calls-nocache"]++;

    // compute score
    // Extracting 'cov_target_target' value
    double cov_target_target = covariance_matrix(target, target);

    double sigma;
    if (parents.empty()) {
        sigma = cov_target_target;
    } else {
        // Building the 'cov_parents_parents' matrix
        const int p_size = parents.size();
        MatrixXd cov_parents_parents(p_size, p_size);
        std::vector<int> parents_vector(parents.begin(), parents.end());
        for (int i = 0; i < p_size; ++i) {
            for (int j = 0; j < p_size; ++j) {
                cov_parents_parents(i, j) =
                        covariance_matrix(parents_vector[i], parents_vector[j]);
            }
        }
        // Building the 'cov_parents_target' vector
        VectorXd cov_parents_target(p_size);
        for (int i = 0; i < p_size; ++i) {
            cov_parents_target(i) = covariance_matrix(parents_vector[i], target);
        }

        VectorXd beta = cov_parents_parents.llt().solve(cov_parents_target);
        sigma = cov_target_target - (cov_parents_target.transpose() * beta).value();
    }
    // Calculating the log-likelihood without the constant
    double log_likelihood_no_constant = -0.5 * n_samples * (1 + std::log(sigma));

    // Calculating the BIC regularization term
    double bic_regularization = 0.5 * std::log(n_samples) * (parents.size() + 1.) * alpha;

    //    double prior_regularization = log_binomial(n_variables - 1, parents.size());
    // makes things worse

    // Calculating the BIC score
    double bic = log_likelihood_no_constant - bic_regularization;
    // bic = bic / n_samples;

    // cache update
    cache_target[parents] = bic;

    return bic;
}
