//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include <set>

class ScorerInterface {
public:
    virtual ~ScorerInterface() = default;

    virtual double local_diff_score(int node, const std::set<int> &parents, int new_parent) const = 0;
};