//
// Created by Achille Nazaret on 11/6/23.
//
#include "XGES.h"
#include "set_ops.h"
#include <iostream>
#include <stack>

using namespace std::chrono;

XGES::XGES(const Eigen::MatrixXd &data, ScorerInterface *scorer)
    : pdag(data.cols()), scorer(scorer), initial_score(scorer->score_pdag(pdag)) {
    n_variables = data.cols();
    n_samples = data.rows();
    total_score = initial_score;
}

XGES::XGES(const Eigen::MatrixXd &data, std::vector<FlatSet> interventions_candidate_variables,
           ScorerInterface *scorer)
    : pdag(data.cols(), interventions_candidate_variables.size()), scorer(scorer),
      initial_score(scorer->score_pdag(pdag)) {
    n_variables = data.cols();
    n_samples = data.rows();
    total_score = initial_score;
    this->interventions_candidate_variables = interventions_candidate_variables;
    n_interventions = interventions_candidate_variables.size();
    variables_candidate_interventions.resize(n_variables);
    for (int i = 0; i < n_variables; ++i) {
        for (int j = 0; j < n_interventions; ++j) {
            if (interventions_candidate_variables[j].find(i) != interventions_candidate_variables[j].end()) {
                variables_candidate_interventions[i].insert(j + n_variables);
            }
        }
    }
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
    auto start_init_inserts = high_resolution_clock::now();
    for (int y = 0; y < n_variables; ++y) {
        // find all possible inserts to y
        find_inserts_to_y(y, candidate_inserts, -1, false);
    }
    statistics["time- init_inserts"] =
            duration_cast<duration<double>>(high_resolution_clock::now() - start_init_inserts).count();
    std::cout << "candidate_inserts.size() = " << candidate_inserts.size() << std::endl;

    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    int duplicate_inserts = 0;
    int n_inserts = 0;
    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    // now do the heuristic
    while (!candidate_inserts.empty() || !candidate_reverses.empty() || !candidate_deletes.empty()) {
        // In order: reverse, delete, insert
        // Apply only one operator per iteration
        edge_modifications.clear();

        if (!candidate_deletes.empty()) {
            // pop the best delete
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();

            // check if it is still valid
            if (pdag.is_delete_valid(best_delete)) {
                // apply the delete
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                // log it
                std::cout << i_operations << ". " << best_delete << std::endl;
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
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                // log it
                std::cout << i_operations << ". " << best_reverse << std::endl;
            } else {
                continue;
            }

        } else if (!candidate_inserts.empty()) {
            // pop the best insert
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();

            n_inserts++;
            if (best_insert.y == last_insert.y && abs(best_insert.score - last_insert.score) < 1e-10 &&
                best_insert.x == last_insert.x && best_insert.T == last_insert.T) {
                duplicate_inserts++;
                continue;
            }
            last_insert = std::move(best_insert);

            // check if it is still valid
            auto start_time = high_resolution_clock::now();
            if (pdag.is_insert_valid(last_insert)) {
                // apply the insert
                statistics["time- is_insert_valid true"] +=
                        duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();
                start_time = high_resolution_clock::now();
                pdag.apply_insert(last_insert, edge_modifications);
                statistics["time- apply_insert"] +=
                        duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();
                total_score += last_insert.score;
                // log it
                std::cout << i_operations << ". " << last_insert << std::endl;
            } else {
                statistics["time- is_insert_valid false"] +=
                        duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();
                continue;
            }
        }
        i_operations++;

        //        candidate_inserts.clear();
        //        candidate_reverses.clear();
        //        candidate_deletes.clear();
        //        for (int node: pdag.get_nodes()) {
        //            find_inserts_to_y(node, candidate_inserts);
        //            find_reverse_to_y(node, candidate_reverses);
        //            find_deletes_to_y(node, candidate_deletes);
        //        }
        //        continue;

        std::set<int> touched_nodes;
        std::set<int> full_insert_to_y;
        for (auto &edge_modification_key_value: edge_modifications) {
            auto &edge = edge_modification_key_value.second;
            touched_nodes.insert(edge.x);
            touched_nodes.insert(edge.y);
            if (edge.is_now_reverse() || edge.is_now_undirected()) {
                full_insert_to_y.insert(edge.x);
                full_insert_to_y.insert(edge.y);
                std::set_intersection(pdag.get_neighbors(edge.x).begin(), pdag.get_neighbors(edge.x).end(),
                                      pdag.get_neighbors(edge.y).begin(), pdag.get_neighbors(edge.y).end(),
                                      std::inserter(full_insert_to_y, full_insert_to_y.begin()));

            } else if (edge.is_now_directed()) {
                if (edge.new_type == EdgeType::DIRECTED_TO_Y) {
                    full_insert_to_y.insert(edge.y);
                } else {
                    full_insert_to_y.insert(edge.x);
                }
                std::set_intersection(pdag.get_neighbors(edge.x).begin(), pdag.get_neighbors(edge.x).end(),
                                      pdag.get_neighbors(edge.y).begin(), pdag.get_neighbors(edge.y).end(),
                                      std::inserter(full_insert_to_y, full_insert_to_y.begin()));
            } else {
                full_insert_to_y.insert(edge.x);
                full_insert_to_y.insert(edge.y);
            }
        }

        auto start_time = high_resolution_clock::now();
        for (auto node: touched_nodes) {
            find_reverse_to_y(node, candidate_reverses);
            find_reverse_from_x(node, candidate_reverses);
            find_deletes_to_y(node, candidate_deletes);
        }
        for (auto node: full_insert_to_y) { find_inserts_to_y(node, candidate_inserts); }

        for (auto &edge_modification: edge_modifications) {
            // symmetric in x and y
            auto &edge = edge_modification.second;
            for (auto target: pdag.get_neighbors(edge.y)) {
                if (full_insert_to_y.find(target) != full_insert_to_y.end()) { continue; }
                find_inserts_to_y(target, candidate_inserts, edge.x);
            }
            for (auto target: pdag.get_neighbors(edge.x)) {
                if (full_insert_to_y.find(target) != full_insert_to_y.end()) { continue; }
                find_inserts_to_y(target, candidate_inserts, edge.y);
            }

            if (edge.is_now_directed()) {
                int x_ = edge.get_source();
                int y_ = edge.get_target();
                if (full_insert_to_y.find(x_) == full_insert_to_y.end()) {
                    for (auto source: pdag.get_adjacent(y_)) {
                        find_inserts_to_y(x_, candidate_inserts, source);
                    }
                }
            }
        }
        statistics["time- update_operators"] +=
                duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();


        //        std::cout << "score=" << total_score << std::endl << std::endl;
    }
    std::cout << "probable_duplicates = " << duplicate_inserts << std::endl;
    std::cout << "n_inserts = " << n_inserts << std::endl;
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
void XGES::find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts, int parent_x,
                             bool low_parent_only, bool positive_only) {
    if (node_is_intervention(y)) { return; }
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y = pdag.get_parents(y);

    std::set<int> possible_parents;

    if (parent_x != -1) {
        possible_parents.insert(parent_x);
    } else {
        // for now: no pre-selection
        auto &nodes = pdag.get_nodes_variables();
        if (n_interventions) {
            // include candidate interventions
            possible_parents.insert(variables_candidate_interventions[y].begin(),
                                    variables_candidate_interventions[y].end());
        }

        // 1. x is not adjacent to y (x ∉ Ad(y))
        std::set_difference(nodes.begin(), nodes.end(), adjacent_y.begin(), adjacent_y.end(),
                            std::inserter(possible_parents, possible_parents.begin()));
        possible_parents.erase(y);// only needed because we don't have pre-selection
    }


    for (int x: possible_parents) {
        if (low_parent_only && x > y) { break; }
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
        FlatSet &effective_parents_y = neighbors_y_adjacent_x;// just renaming it, no copy necessary
        effective_parents_y.insert(parents_y.begin(), parents_y.end());
        // <T: set of nodes, iterator over neighbors_y_not_adjacent_x, effective_parents: set of nodes>
        std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
        // we know that T = {} is valid
        stack.emplace(FlatSet{}, neighbors_y_not_adjacent_x.begin(), effective_parents_y);

        while (!stack.empty()) {
            auto top = std::move(stack.top());
            stack.pop();
            auto &T = std::get<0>(top);
            auto it = std::get<1>(top);
            auto &effective_parents = std::get<2>(top);

            // change if we parallelize
            auto start = high_resolution_clock::now();
            double score = scorer->score_insert(y, effective_parents, x);
            statistics["time- score_insert"] +=
                    duration_cast<duration<double>>(high_resolution_clock::now() - start).count();
            if (score > 0 || !positive_only) {
                // using move(T)/move(effective_parents) should also work even though we look them up
                // later. but we play it safe for now.
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
                    auto T_prime = T;
                    T_prime.insert(z);
                    auto effective_parents_prime = effective_parents;
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
        std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
        // we know that O = {} is valid
        FlatSet effective_parents_init;
        effective_parents_init.reserve(parents_y.size() + neighbors_y_adjacent_x.size() + 1);
        // note: Chickering Corollary 18 is incorrect. Pa(y) might not contain x, it has to be added.
        union_with_single_element(parents_y, x, effective_parents_init);
        stack.emplace(FlatSet{}, neighbors_y_adjacent_x.begin(), effective_parents_init);

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
                    auto O_prime = O;
                    O_prime.insert(z);
                    auto effective_parents_prime = effective_parents;
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
        auto &parents_x = pdag.get_parents(x);
        std::vector<Insert> candidate_inserts;
        find_inserts_to_y(y, candidate_inserts, x, false, false);

        for (auto insert: candidate_inserts) {
            // change if we parallelize
            double score = insert.score + scorer->score_delete(x, parents_x, y);

            if (score > 0) {
                candidate_reverses.emplace_back(insert, score, parents_x);
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
        find_inserts_to_y(y, candidate_inserts, x, false, false);

        for (auto insert: candidate_inserts) {
            // change if we parallelize
            double score = insert.score + scorer->score_delete(x, parents_x, y);

            if (score > 0) {
                candidate_reverses.emplace_back(insert, score, parents_x);
                std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
            }
        }
    }
}


const PDAG &XGES::get_pdag() const { return pdag; }

double XGES::get_score() const { return total_score; }
