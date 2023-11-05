//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include <set>
#include "ScorerInterface.h"

class Insert {
public:
    Insert(int x, int y, const std::set<int> &T);

    double compute_score(const ScorerInterface &scorer);

    friend std::ostream &operator<<(std::ostream &os, const Insert &obj);

    friend class PDAG;

private:
    int x, y;
    std::set<int> T;
    std::set<int> effective_parents; // = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    double score = 0;
};


