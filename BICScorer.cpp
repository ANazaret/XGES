//
// Created by Achille Nazaret on 11/7/23.
//
#include "BICScorer.h"
#include <iostream>

MatrixXd compute_covariance(const MatrixXd &data) {
    int n_variables = data.cols();
    int n_samples = data.rows();

    MatrixXd centered = data.rowwise() - data.colwise().mean();
    MatrixXd covariance_matrix = (centered.adjoint() * centered);
    return covariance_matrix / n_samples;
}

MatrixXd compute_covariance(const MatrixXd &data, const VectorXi &interventions_index) {
    int n_variables = data.cols();
    int n_interventions = interventions_index.maxCoeff() + 1;
    int n_samples = data.rows();
    // With interventions, the covariance matrix is of size (n_variables + n_interventions)^2
    // Top left is the covariance matrix of the variables
    // Top right and bottom left are the covariances between variables and interventions
    // Bottom right is the covariance matrix of the interventions

    MatrixXd covariance_matrix(n_variables + n_interventions, n_variables + n_interventions);

    MatrixXd centered = data.rowwise() - data.colwise().mean();
    covariance_matrix.topLeftCorner(n_variables, n_variables) = (centered.adjoint() * centered) / n_samples;


    MatrixXd cov_obs_inter = MatrixXd::Zero(n_interventions, n_variables);
    VectorXd means_inter = VectorXd::Zero(n_interventions);
    MatrixXd cov_inter_inter = MatrixXd::Zero(n_interventions, n_interventions);

    for (int i = 0; i < n_samples; ++i) {
        int i_intervention = interventions_index(i);
        if (i_intervention == -1) { continue; }
        cov_obs_inter.row(i_intervention) += centered.row(i);
        means_inter(i_intervention) += 1;
        cov_inter_inter(i_intervention, i_intervention) += 1;
    }
    cov_obs_inter /= n_samples;
    means_inter /= n_samples;
    cov_inter_inter /= n_samples;
    cov_inter_inter -= means_inter * means_inter.transpose();// not sure this works

    covariance_matrix.topRightCorner(n_variables, n_interventions) = cov_obs_inter.transpose();
    covariance_matrix.bottomLeftCorner(n_interventions, n_variables) = cov_obs_inter;
    covariance_matrix.bottomRightCorner(n_interventions, n_interventions) = cov_inter_inter;

    return covariance_matrix;
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

double log_binomial(int n, int k) {
    // use the fact that log(n choose k) = log(n!) - log(k!) - log((n-k)!)
    // and log-gamma(n+1) = log(n!)
    double log_n_choose_k = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
    return log_n_choose_k;
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

    //    double prior_regularization = log_binomial(n_variables - 1, parents.size());
    // makes things worse

    // Calculating the BIC score
    double bic = log_likelihood_no_constant - bic_regularization;

    // cache update
    cache_target[parents] = bic;

    return bic;
}