//
// Created by Achille Nazaret on 11/7/23.
//
#include "BICScorer.h"
#include <iostream>

MatrixXd compute_covariance(const MatrixXd &data) {
    int n_variables = data.cols();
    int n_samples = data.rows();
    MatrixXd covariance_matrix(n_variables, n_variables);
    for (int i = 0; i < n_variables; ++i) {
        for (int j = 0; j < n_variables; ++j) {
            double covariance = (data.col(i).array() - data.col(i).mean()).matrix().transpose() *
                                (data.col(j).array() - data.col(j).mean()).matrix();
            covariance_matrix(i, j) = covariance;
        }
    }
    // todo: benchmark, benchmark, benchmark (move n_samples in the loop, and ultimately try
    // vectorization)
    return covariance_matrix / n_samples;
}

BICScorer::BICScorer(const Eigen::MatrixXd &data, double alpha)
    : data(data), alpha(alpha), covariance_matrix(compute_covariance(data)) {
    n_variables = data.cols();
    n_samples = data.rows();
}

double BICScorer::local_diff_score(int target, const std::set<int> &parents, int new_parent) const {
    //compute local_score(target, parents U {new_parent}) - local_score(target, parents)
    // can be optimized to not create a new set
    std::set<int> parents_with_new_parent = parents;
    parents_with_new_parent.insert(new_parent);
    double score_with_new_parent = local_score(target, parents_with_new_parent);
    double score_without_new_parent = local_score(target, parents);
    return score_with_new_parent - score_without_new_parent;
}

double BICScorer::local_score(int target, const std::set<int> &parents) const {
    // cache lookup
    auto it = cache.find(std::make_pair(target, parents));
    if (it != cache.end()) { return it->second; }
    // compute score
    // Extracting 'cov_target_target' value
    double cov_target_target = covariance_matrix(target, target);

    double sigma;
    if (parents.empty()) {
        sigma = cov_target_target;
    } else {
        // Building the 'cov_parents_parents' matrix
        int p_size = parents.size();
        Eigen::MatrixXd cov_parents_parents(p_size, p_size);
        std::vector<int> parents_vector(parents.begin(), parents.end());
        for (int i = 0; i < p_size; ++i) {
            for (int j = 0; j < p_size; ++j) {
                cov_parents_parents(i, j) = covariance_matrix(parents_vector[i], parents_vector[j]);
            }
        }
        // Building the 'cov_parents_target' vector
        Eigen::VectorXd cov_parents_target(p_size);
        for (int i = 0; i < p_size; ++i) {
            cov_parents_target(i) = covariance_matrix(parents_vector[i], target);
        }

        Eigen::VectorXd beta = cov_parents_parents.llt().solve(cov_parents_target);
        sigma = cov_target_target - (cov_parents_target.transpose() * beta).value();
    }
    // Calculating the log-likelihood without the constant
    double log_likelihood_no_constant = -0.5 * n_samples * (1 + std::log(sigma));

    // Calculating the BIC regularization term
    double bic_regularization = 0.5 * std::log(n_samples) * (parents.size() + 1) * alpha;

    // Calculating the BIC score
    double bic = log_likelihood_no_constant - bic_regularization;

    return bic;
}