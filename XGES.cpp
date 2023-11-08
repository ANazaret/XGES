//
// Created by Achille Nazaret on 11/6/23.
//
#include "XGES.h"
#include "set_ops.h"
#include <iostream>
#include <stack>

XGES::XGES(const Eigen::MatrixXd &data, ScorerInterface *scorer) : pdag(data.cols()), scorer(scorer) {
    n_variables = data.cols();
    n_samples = data.rows();
}


void XGES::fit_heuristic() {
    // pre-selection

    // heuristic turn
    heuristic_turn_delete_insert();
}

void XGES::heuristic_turn_delete_insert() {
    std::vector<Insert> candidate_inserts;
    candidate_inserts.reserve(100 * n_variables);

    // init the candidate inserts
    std::cout << "n_variables = " << n_variables << std::endl;
    for (int y = 0; y < n_variables; ++y) {
        // find all possible inserts to y
        find_inserts_to_y(y, candidate_inserts);
    }
    std::cout << "candidate_inserts.size() = " << candidate_inserts.size() << std::endl;


    // now do the heuristic
    while (!candidate_inserts.empty()) {
        // pop the best insert
        std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
        auto best_insert = std::move(candidate_inserts.back());
        candidate_inserts.pop_back();

        // check if it is still valid
        if (pdag.is_insert_valid(best_insert)) {
            // apply the insert
            pdag.apply_insert(best_insert);
            total_score += best_insert.score;
            // log it
            std::cout << "Insert: " << best_insert << std::endl;
            // update the candidate inserts
            // 1. remove all inserts to y
            //            candidate_inserts.erase(
            //                    std::remove_if(candidate_inserts.begin(), candidate_inserts.end(),
            //                                   [best_insert](const Insert &insert) { return insert.y == best_insert.y; }),
            //                    candidate_inserts.end());
            // 2. find all possible inserts to y
            find_inserts_to_y(best_insert.y, candidate_inserts);
            find_inserts_to_y(best_insert.x, candidate_inserts);

            std::cout << pdag << std::endl;
        }
    }
}


/**
 * Find all possible inserts to y.
 *
 * The candidate inserts (x, y, T) are such that:
 *  1. x is not adjacent to y (x ∉ Ad(y))
 *  2. T ⊆ Ne(y) \ Ad(x)
 *  3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 *
 * @param y
 * @param candidate_inserts
 */
void XGES::find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts) {
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y = pdag.get_parents(y);

    // for now: no pre-selection
    auto &nodes = pdag.get_nodes();
    std::set<int> possible_parents;
    // 1. x is not adjacent to y (x ∉ Ad(y))
    std::set_difference(nodes.begin(), nodes.end(), adjacent_y.begin(), adjacent_y.end(),
                        std::inserter(possible_parents, possible_parents.begin()));
    possible_parents.erase(y);// only needed because we don't have pre-selection

    for (int x: possible_parents) {
        if (x == 1 && y == 4) {
            std::cout << "x = " << x;
            std::cout << ", y = " << y << std::endl;
        }
        // 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
        // So in particular, [Ne(y) ∩ Ad(x)] must be a clique.
        auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
        if (!pdag.is_clique(neighbors_y_adjacent_x)) { continue; }

        // 2. T ⊆ Ne(y) \ Ad(x)
        auto neighbors_y_not_adjacent_x = pdag.get_neighbors_not_adjacent(y, x);

        // We enumerate all T ⊆ Ne(y) \ Ad(x) such that [Ne(y) ∩ Ad(x)] ∪ T is a clique (noted C(x,y,T)).
        // If the condition fails for some T, it will fail for all its supersets.
        // Hence, we enumerate the T in inclusion order, to minimize the number of T we need to check.
        // We simulate a recursive search using a stack.
        // The stack contains valid T, and an iterator for the next entries in neighbors_y_not_adjacent_x to consider.

        // The effective parents_y are [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y).
        std::set<int> effective_parents_y;
        std::set_union(neighbors_y_adjacent_x.begin(), neighbors_y_adjacent_x.end(), parents_y.begin(),
                       parents_y.end(), std::inserter(effective_parents_y, effective_parents_y.begin()));
        // <set of nodes, iterator over neighbors_y_not_adjacent_x, set of effective_parents>
        std::stack<std::tuple<std::set<int>, std::set<int>::iterator, std::set<int>>> stack;
        // we know that T = {} is valid
        stack.emplace(std::set<int>{}, neighbors_y_not_adjacent_x.begin(), effective_parents_y);

        while (!stack.empty()) {
            auto top = std::move(stack.top());
            stack.pop();
            auto T = std::get<0>(top);
            auto it = std::get<1>(top);
            auto effective_parents = std::get<2>(top);

            // change if we parallelize
            double score = scorer->score_insert(y, effective_parents, x);
            if (score > 0) {
                candidate_inserts.emplace_back(x, y, T, score, effective_parents);
                std::push_heap(candidate_inserts.begin(), candidate_inserts.end());
            }

            // Look for other candidate T using the iterator, which gives us the next elements z to consider.
            while (it != neighbors_y_not_adjacent_x.end()) {
                // We define T' = T ∪ {z} and we check C(x,y,T') is a clique.
                // Since C(x,y,T) was a clique, we only need to check that z is adjacent to all nodes in C(x,y,T).
                auto z = *it;
                ++it;
                auto &adjacent_z = pdag.get_adjacent(z);
                // We check that C(x,y,T) ⊆ Ad(z); i.e. T ⊆ Ad(z) and [Ne(y) ∩ Ad(x)] ⊆ Ad(z).
                if (std::includes(adjacent_z.begin(), adjacent_z.end(), T.begin(), T.end()) &&
                    std::includes(adjacent_z.begin(), adjacent_z.end(), neighbors_y_adjacent_x.begin(),
                                  neighbors_y_adjacent_x.end())) {
                    // T' is a candidate
                    std::set<int> T_prime = T;
                    T_prime.insert(z);
                    std::set<int> effective_parents_prime = effective_parents;
                    effective_parents_prime.insert(z);
                    stack.emplace(std::move(T_prime), it, std::move(effective_parents_prime));
                }
            }
        }
    }
}
