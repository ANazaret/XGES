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

    std::vector<Reverse> candidate_reverses;
    candidate_inserts.reserve(n_variables);

    std::vector<Delete> candidate_deletes;
    candidate_deletes.reserve(n_variables);


    // init the candidate inserts
    std::cout << "n_variables = " << n_variables << std::endl;
    clock_t start_init_inserts = clock();
    for (int y = 0; y < n_variables; ++y) {
        // find all possible inserts to y
        find_inserts_to_y(y, candidate_inserts, -1, true);
    }
    statistics["time- init_inserts"] = (double) (clock() - start_init_inserts) / CLOCKS_PER_SEC;
    std::cout << "candidate_inserts.size() = " << candidate_inserts.size() << std::endl;

    int i_operations = 1;
    // now do the heuristic
    while (!candidate_inserts.empty() || !candidate_reverses.empty() || !candidate_deletes.empty()) {
        // In order: reverse, delete, insert
        // Apply only one operator per iteration
        int x = -1;
        int y = -1;

        std::set<Edge> changed_edges;

        if (!candidate_deletes.empty()) {
            // pop the best delete
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();

            // check if it is still valid
            if (pdag.is_delete_valid(best_delete)) {
                // apply the delete
                pdag.apply_delete(best_delete, changed_edges);
                total_score += best_delete.score;
                // log it
                std::cout << i_operations << ". " << best_delete << std::endl;

                x = best_delete.x;
                y = best_delete.y;
            } else {
                continue;
            }

        } else if (!candidate_reverses.empty()) {
            // pop the best reverse
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();

            // check if it is still valid
            if (pdag.is_reverse_valid(best_reverse)) {
                // apply the reverse
                pdag.apply_reverse(best_reverse, changed_edges);
                total_score += best_reverse.score;
                // log it
                std::cout << i_operations << ". " << best_reverse << std::endl;
                x = best_reverse.insert.x;
                y = best_reverse.insert.y;
            } else {
                continue;
            }

        } else if (!candidate_inserts.empty()) {
            // pop the best insert
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();

            // check if it is still valid
            if (pdag.is_insert_valid(best_insert)) {
                // apply the insert
                pdag.apply_insert(best_insert, changed_edges);
                total_score += best_insert.score;
                // log it
                std::cout << i_operations << ". " << best_insert << std::endl;

                x = best_insert.x;
                y = best_insert.y;
            } else {
                continue;
            }
        }


        // Temporary logic to update operators
        std::set<int> touched_nodes;
        for (auto &edge: changed_edges) {
            touched_nodes.insert(edge.x);
            touched_nodes.insert(edge.y);
        }

        for (auto node: touched_nodes) {
            find_inserts_to_y(node, candidate_inserts);
            find_reverse_to_y(node, candidate_reverses);
            find_reverse_from_x(node, candidate_reverses);
            find_deletes_to_y(node, candidate_deletes);
        }

        for (auto target: pdag.get_neighbors(y)) {
            if (touched_nodes.find(target) != touched_nodes.end()) { continue; }
            find_inserts_to_y(target, candidate_inserts, x);
        }

        for (auto target: pdag.get_neighbors(x)) {
            if (touched_nodes.find(target) != touched_nodes.end()) { continue; }
            find_inserts_to_y(target, candidate_inserts, y);
        }

        i_operations++;

        std::cout << "score=" << total_score << std::endl << std::endl;
    }
}


/**
 * Find all possible inserts to y.
 *
 * The candidate inserts (x, y, T) are such that:
 *  1. x is not adjacent to y (x ∉ Ad(y))
 *  2. T ⊆ Ne(y) \ Ad(x)
 *  3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 *  Not enforced at that stage: [Ne(y) ∩ Ad(x)] ∪ T blocks all semi-directed paths from y to x
 *
 * @param y
 * @param candidate_inserts
 */
void XGES::find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts, int parent_x) {
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y = pdag.get_parents(y);

    std::set<int> possible_parents;

    if (parent_x != -1) {
        possible_parents.insert(parent_x);
    } else {
        // for now: no pre-selection
        auto &nodes = pdag.get_nodes();
        // 1. x is not adjacent to y (x ∉ Ad(y))
        std::set_difference(nodes.begin(), nodes.end(), adjacent_y.begin(), adjacent_y.end(),
                            std::inserter(possible_parents, possible_parents.begin()));
        possible_parents.erase(y);// only needed because we don't have pre-selection
    }


    for (int x: possible_parents) {
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
        // <T: set of nodes, iterator over neighbors_y_not_adjacent_x, effective_parents: set of nodes>
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
            clock_t start = clock();
            double score = scorer->score_insert(y, effective_parents, x);
            statistics["time- score_insert"] += (double) (clock() - start) / CLOCKS_PER_SEC;
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

/** Find all possible deletes to y.
 *
 * The candidate deletes (x, y, H) are such that:
 *  1. x is a parent or a neighbor of y
 *  2. H ⊆ [Ne(y) ∩ Ad(x)]  (alt. O ⊆ [Ne(y) ∩ Ad(x)] -- the complement of H in [Ne(y) ∩ Ad(x)])
 *  3. [Ne(y) ∩ Ad(x)] \ H is a clique (alt. O is a clique)
 *
 *  The effective parents_y are Pa(y) ∪ [Ne(y) ∩ Ad(x)] \ H (alt. Pa(y) ∪ O)
 *
 * @param y
 * @param candidate_deletes
 */
void XGES::find_deletes_to_y(int y, std::vector<Delete> &candidate_deletes) {
    auto &neighbors_y = pdag.get_neighbors(y);
    auto &parents_y = pdag.get_parents(y);

    std::vector<int> possible_x;
    possible_x.insert(possible_x.end(), parents_y.begin(), parents_y.end());
    possible_x.insert(possible_x.end(), neighbors_y.begin(), neighbors_y.end());

    for (int x: possible_x) {
        auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);

        // find all possible O ⊆ [Ne(y) ∩ Ad(x)] that are cliques, quite similar to the inserts but simpler
        // <O set of nodes, iterator over neighbors_y_adjacent_x, set of effective_parents>
        std::stack<std::tuple<std::set<int>, std::set<int>::iterator, std::set<int>>> stack;
        // we know that O = {} is valid
        stack.emplace(std::set<int>{}, neighbors_y_adjacent_x.begin(), parents_y);

        while (!stack.empty()) {
            auto top = std::move(stack.top());
            stack.pop();
            auto O = std::get<0>(top);
            auto it = std::get<1>(top);
            auto effective_parents = std::get<2>(top);

            // change if we parallelize
            double score = scorer->score_delete(y, effective_parents, x);
            if (score > 0) {
                candidate_deletes.emplace_back(x, y, O, score, effective_parents);
                std::push_heap(candidate_deletes.begin(), candidate_deletes.end());
            }

            // Look for other candidate O using the iterator, which gives us the next elements z to consider.
            while (it != neighbors_y_adjacent_x.end()) {
                // We define O' = O ∪ {z} and we check if O' is a clique.
                // Since O was a clique, we only need to check that z is adjacent to all nodes in O.
                auto z = *it;
                ++it;
                auto &adjacent_z = pdag.get_adjacent(z);
                // We check that O ⊆ Ad(z)
                if (std::includes(adjacent_z.begin(), adjacent_z.end(), O.begin(), O.end())) {
                    // O' is a candidate
                    std::set<int> O_prime = O;
                    O_prime.insert(z);
                    std::set<int> effective_parents_prime = effective_parents;
                    effective_parents_prime.insert(z);
                    stack.emplace(std::move(O_prime), it, std::move(effective_parents_prime));
                }
            }
        }
    }
}

/**
 * Find if I can turn x ← y into x → y.
 * The candidate turn (x, y, T) are such that:
 *  2. T ⊆ Ne(y) \ Ad(x)
 *  3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 *  Not enforced at that stage: [Ne(y) ∩ Ad(x)] ∪ T blocks all semi-directed paths from y to x (≠ than x ← y)
 *
 * @param y
 * @param candidate_inserts
 */

// can i just do insert x -> y?
void XGES::find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses) {
    // look for all possible x ← y
    auto &children_y = pdag.get_children(y);

    for (int x: children_y) {

        std::vector<Insert> candidate_inserts;
        find_inserts_to_y(y, candidate_inserts, x);

        for (auto insert: candidate_inserts) {
            // change if we parallelize
            double score = insert.score + scorer->score_delete(x, pdag.get_parents(x), y);

            if (score > 0) {
                candidate_reverses.emplace_back(insert, score);
                std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
            }
        }
    }
}

void XGES::find_reverse_from_x(int x, std::vector<Reverse> &candidate_reverses) {
    // look for all possible x ← y
    auto &parents_x = pdag.get_parents(x);


    for (int y: parents_x) {

        std::vector<Insert> candidate_inserts;
        find_inserts_to_y(y, candidate_inserts, x);

        for (auto insert: candidate_inserts) {
            // change if we parallelize
            double score = insert.score + scorer->score_delete(x, pdag.get_parents(x), y);

            if (score > 0) {
                candidate_reverses.emplace_back(insert, score);
                std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
            }
        }
    }
}


const PDAG &XGES::get_pdag() const { return pdag; }

double XGES::get_score() const { return total_score; }