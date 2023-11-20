//
// Created by Achille Nazaret on 11/7/23.
//

#pragma once

#include "ScorerInterface.h"
#include <Eigen/Dense>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using Eigen::MatrixXd;


// Specialize std::hash for std::set<int>
namespace std {
    template<>
    struct hash<set<int>> {
        size_t operator()(const set<int> &s) const {
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
    const MatrixXd &data;
    double alpha;
    const MatrixXd covariance_matrix;

    // implement a cache for the local_score function: map from (target, parents) to score
    //    std::map<std::pair<int, std::set<int>>, double> cache;
    std::vector<std::unordered_map<std::set<int>, double>> cache;


public:
    BICScorer(const MatrixXd &data, double alpha);
    double local_score(int target, const std::set<int> &parents) override;
};
