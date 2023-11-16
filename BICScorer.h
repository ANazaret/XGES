//
// Created by Achille Nazaret on 11/7/23.
//

#pragma once

#include "ScorerInterface.h"
#include <Eigen/Dense>
#include <map>

using Eigen::MatrixXd;

class BICScorer : public ScorerInterface {
private:
    int n_variables;
    int n_samples;
    const MatrixXd &data;
    double alpha;
    const MatrixXd covariance_matrix;

    // implement a cache for the local_score function: map from (target, parents) to score
    std::map<std::pair<int, std::set<int>>, double> cache;


public:
    BICScorer(const MatrixXd &data, double alpha);
    double local_score(int target, const std::set<int> &parents) const override;
};
