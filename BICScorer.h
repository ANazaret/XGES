//
// Created by Achille Nazaret on 11/7/23.
//

#pragma once

#include "ScorerInterface.h"
#define EIGEN_USE_BLAS
#include "Eigen/Dense"
#include <map>
#include <unordered_map>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// Specialize std::hash for std::set<int>
namespace std {
    template<>
    struct hash<FlatSet> {
        size_t operator()(const FlatSet &s) const {
            size_t hash_value = 0;
            for (const int &elem: s) {
                hash_value ^= hash<int>{}(elem) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
            }
            return hash_value;
        }
    };
}// namespace std

class BICScorer : public ScorerInterface {
private:
    int n_variables;
    int n_samples;
    int n_interventions = 0;

    const MatrixXd &data;
    const VectorXi interventions_index = VectorXi::Zero(0);
    double alpha;
    const MatrixXd covariance_matrix;

    // implement a cache for the local_score function: map from (target, parents) to score
    //    std::map<std::pair<int, std::set<int>>, double> cache;
    std::vector<std::unordered_map<FlatSet, double>> cache;


public:
    BICScorer(const MatrixXd &data, double alpha);
    BICScorer(const MatrixXd &data, const VectorXi &interventions_index, double alpha);
    double local_score(int target, const FlatSet &parents) override;
};
