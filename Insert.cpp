//
// Created by Achille Nazaret on 11/4/23.
//

#include "Insert.h"

Insert::Insert(int x, int y, const std::set<int> &T) : x(x), y(y), T(T) {
    // need to compute effective parents

}

double Insert::compute_score(const ScorerInterface &scorer) {
    score = scorer.local_diff_score(y, T, x);
    return score;
}
