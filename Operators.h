//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include "ScorerInterface.h"
#include <iostream>
#include <set>

class Insert {
public:
    Insert(int x, int y, const std::set<int> &T, double score, std::set<int> effective_parents);

    friend std::ostream &operator<<(std::ostream &os, const Insert &obj);

    friend class PDAG;

    friend class XGES;

    bool operator<(const Insert &rhs) const { return score < rhs.score; }


private:
    int x, y;
    std::set<int> T;
    std::set<int> effective_parents;// = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    double score = 0;
};


class Delete {
public:
    Delete(int x, int y, const std::set<int> &O, double score, std::set<int> effective_parents);

    friend std::ostream &operator<<(std::ostream &os, const Delete &obj);

    friend class PDAG;

    friend class XGES;

    bool operator<(const Delete &rhs) const { return score < rhs.score; }


private:
    int x, y;
    std::set<int> H;
    // std::set<int> effective_parents;// = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    std::set<int> O;
    //std::set<int> effective_parents; // =
    double score = 0;
};
    double score = 0;
};