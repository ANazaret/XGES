//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once


#include "set_ops.h"
#include <iostream>
#include <set>

class Reverse;


class Insert {
public:
    Insert(int x, int y, const FlatSet &T, double score, const FlatSet &effective_parents);

    friend std::ostream &operator<<(std::ostream &os, const Insert &obj);
    friend std::ostream &operator<<(std::ostream &os, const Reverse &obj);

    friend class PDAG;

    friend class XGES;

    bool operator<(const Insert &rhs) const { return score < rhs.score; }


private:
    int x, y;
    FlatSet T;
    FlatSet effective_parents;// = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    double score = 0;
};

class Delete {
public:
    Delete(int x, int y, const FlatSet &O, double score, const FlatSet &effective_parents, bool directed);

    friend std::ostream &operator<<(std::ostream &os, const Delete &obj);

    friend class PDAG;

    friend class XGES;

    bool operator<(const Delete &rhs) const { return score < rhs.score; }


private:
    int x, y;
    bool directed;
    FlatSet O;
    FlatSet effective_parents;// O ∪ Pa(y)
    double score = 0;
};


class Reverse {
public:
    Reverse(int x, int y, const FlatSet &T, double score, const FlatSet &effective_parents,
            const FlatSet &parents_x);
    Reverse(Insert insert, double score, const FlatSet &parents_x);

    friend std::ostream &operator<<(std::ostream &os, const Reverse &obj);

    friend class PDAG;

    friend class XGES;

    bool operator<(const Reverse &rhs) const { return score < rhs.score; }


private:
    Insert insert;
    double score = 0;
    FlatSet parents_x;
};