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

MatrixXd compute_covariance(const MatrixXd &data, const VectorXi &interventions_index) {
    int n_variables = data.cols();
    int n_interventions = interventions_index.maxCoeff() + 1;
    int n_samples = data.rows();
    Eigen::VectorXd means = data.colwise().mean();

    MatrixXd covariance_matrix(n_variables + n_interventions, n_variables + n_interventions);
    for (int i = 0; i < n_variables + n_interventions; ++i) {
        for (int j = 0; j <= i; ++j) {
            Eigen::VectorXd x;
            Eigen::VectorXd y;
            if (i < n_variables) {
                x = data.col(i).array() - means(i);
            } else {
                int i_intervention = i - n_variables;
                x = Eigen::VectorXd::Zero(n_samples);
                // set x to 1 for samples where the intervention is active
                for (int k = 0; k < n_samples; ++k) {
                    if (interventions_index(k) == i_intervention) { x(k) = 1; }
                }
                x.array() -= x.mean();
            }
            if (j < n_variables) {
                y = data.col(j).array() - means(j);
            } else {
                int j_intervention = j - n_variables;
                y = Eigen::VectorXd::Zero(n_samples);
                // set y to 1 for samples where the intervention is active
                for (int k = 0; k < n_samples; ++k) {
                    if (interventions_index(k) == j_intervention) { y(k) = 1; }
                }
                y.array() -= y.mean();
            }
            double covariance = x.matrix().transpose() * y.matrix();
            covariance_matrix(i, j) = covariance;
            covariance_matrix(j, i) = covariance;
        }
    }
    // todo: benchmark, benchmark, benchmark, try vectorization
    return covariance_matrix / n_samples;
}

BICScorer::BICScorer(const Eigen::MatrixXd &data, double alpha)
    : data(data), alpha(alpha), covariance_matrix(compute_covariance(data)) {
    n_variables = data.cols();
    n_samples = data.rows();
    cache.resize(n_variables);
}

BICScorer::BICScorer(const Eigen::MatrixXd &data, const Eigen::VectorXi &interventions_index, double alpha)
    : data(data), alpha(alpha), covariance_matrix(compute_covariance(data, interventions_index)),
      interventions_index(interventions_index) {
    n_variables = data.cols();
    n_interventions = interventions_index.maxCoeff() + 1;
    n_samples = data.rows();
    cache.resize(n_variables);
}

double BICScorer::local_score(int target, const FlatSet &parents) {
    // cache lookup
    auto &cache_target = cache[target];
    auto it = cache_target.find(parents);
    if (it != cache_target.end()) { return it->second; }

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

    // cache update
    cache_target[parents] = bic;

    return bic;
}