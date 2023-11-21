//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include "Operators.h"
#include <set>

class ScorerInterface {
public:
    virtual ~ScorerInterface() = default;

    virtual double local_score(int target, const FlatSet &parents) = 0;

    double score_insert(int target, const FlatSet &parents, int parent_to_add) {
        double score_without_new_parent = local_score(target, parents);
        FlatSet parents_with_new_parent;
        parents_with_new_parent.reserve(parents.size() + 1);
        union_with_single_element(parents, parent_to_add, parents_with_new_parent);

        double score_with_new_parent = local_score(target, parents_with_new_parent);
        return score_with_new_parent - score_without_new_parent;
    }

    double score_delete(int target, const FlatSet &parents, int parent_to_remove) {
        double score_with_old_parent = local_score(target, parents);
        FlatSet parents_without_old_parent;
        parents_without_old_parent.reserve(parents.size() - 1);
        for (auto p: parents) {
            if (p != parent_to_remove) { parents_without_old_parent.insert(p); }
        }
        double score_without_old_parent = local_score(target, parents_without_old_parent);
        return score_without_old_parent - score_with_old_parent;
    }
};