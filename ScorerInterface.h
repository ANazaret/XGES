//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include <set>

class ScorerInterface {
public:
    virtual ~ScorerInterface() = default;

    virtual double local_score(int target, const std::set<int> &parents) const = 0;

    double score_insert(int target, const std::set<int> &parents, int parent_to_add) {
        std::set<int> parents_with_new_parent = parents;
        parents_with_new_parent.insert(parent_to_add);
        double score_with_new_parent = local_score(target, parents_with_new_parent);
        double score_without_new_parent = local_score(target, parents);
        return score_with_new_parent - score_without_new_parent;
    }

    double score_delete(int target, const std::set<int> &parents, int parent_to_remove) {
        std::set<int> parents_without_old_parent = parents;
        parents_without_old_parent.erase(parent_to_remove);
        double score_with_old_parent = local_score(target, parents);
        double score_without_old_parent = local_score(target, parents_without_old_parent);
        return score_without_old_parent - score_with_old_parent;
    }
};