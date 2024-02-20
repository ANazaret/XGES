//
// Created by Achille Nazaret on 11/6/23.
//
#include "XGES.h"
#include "set_ops.h"
#include "utils.h"
#include "EdgeQueueSet.h"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

#include <stack>

using namespace std::chrono;

XGES::XGES(const int n_variables, ScorerInterface *scorer)
    : n_variables(n_variables), scorer(scorer), pdag(n_variables),
      initial_score(scorer->score_pdag(pdag)) {
    total_score = initial_score;
    _logger = spdlog::get("stdout_logger");
}

XGES::XGES(const XGES &other)
    : n_variables(other.n_variables), scorer(other.scorer), pdag(other.pdag),
      initial_score(other.initial_score), total_score(other.total_score),
      _logger(other._logger) {
    // The pointer to the ground truth PDAG is not copied. Change if needed.
}


void XGES::fit_xges(bool extended_search) {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_inserts.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    heuristic_xges0(candidate_inserts, candidate_reverses, candidate_deletes,
                    unblocked_paths_map, true);

    if (extended_search) { block_each_edge_and_research(unblocked_paths_map); }
}


void XGES::block_each_edge_and_research(UnblockedPathsMap &unblocked_paths_map) {
    std::vector<Delete> all_edge_deletes;
    bool deletes_of_pdag_are_updated = false;
    int index = 0;

    while (index < all_edge_deletes.size() || !deletes_of_pdag_are_updated) {
        if (index >= all_edge_deletes.size()) {
            all_edge_deletes.clear();
            for (const int y: pdag.get_nodes_variables()) {
                find_deletes_to_y(y, all_edge_deletes, false);
            }
            deletes_of_pdag_are_updated = true;
            index = 0;
        }

        auto delete_ = all_edge_deletes[index++];
        if (!pdag.is_delete_valid(delete_)) { continue; }

        // Apply the delete
        XGES xges_copy(*this);
        EdgeModificationsMap edge_modifications;
        std::vector<Insert> candidate_inserts;
        std::vector<Reverse> candidate_reverses;
        std::vector<Delete> candidate_deletes;
        xges_copy.pdag.apply_delete(delete_, edge_modifications);
        xges_copy.total_score += delete_.score;
        xges_copy.pdag.add_forbidden_insert(delete_.x, delete_.y);
        UnblockedPathsMap blocked_paths_map_copy = unblocked_paths_map;
        _logger->debug("EXTENDED SEARCH: {}", delete_);
        for (auto &[fst, snd]: edge_modifications) { _logger->trace("\tEdge {}", snd); }


        // xges_copy.update_operator_candidates(edge_modifications, candidate_inserts,
        //                                      candidate_reverses, candidate_deletes);

        xges_copy.update_operator_candidates_v2(edge_modifications, candidate_inserts,
                                                candidate_reverses, candidate_deletes,
                                                blocked_paths_map_copy);

        // xges_copy.update_operator_candidates_naive(candidate_inserts, candidate_reverses,
        //                                            candidate_deletes);

        xges_copy.heuristic_xges0(candidate_inserts, candidate_reverses,
                                  candidate_deletes, blocked_paths_map_copy, false);
        if (pdag == xges_copy.pdag) { continue; }
        if (xges_copy.total_score - total_score > 1e-7) {
            _logger->debug("EXTENDED SEARCH ACCEPTED: {} {}", delete_,
                           xges_copy.total_score);
            total_score = xges_copy.total_score;
            pdag = xges_copy.pdag;
            unblocked_paths_map = std::move(blocked_paths_map_copy);
            deletes_of_pdag_are_updated = false;
        } else {
            _logger->debug("EXTENDED SEARCH REJECTED: {} {}", delete_,
                           xges_copy.total_score);
        }
    }
}


void XGES::heuristic_xges0(std::vector<Insert> &candidate_inserts,
                           std::vector<Reverse> &candidate_reverses,
                           std::vector<Delete> &candidate_deletes,
                           UnblockedPathsMap &unblocked_paths_map,
                           bool initialize_inserts) {

    if (initialize_inserts) {
        // find all possible inserts
        auto start_init_inserts = high_resolution_clock::now();
        for (int y = 0; y < n_variables; ++y) {
            find_inserts_to_y(y, candidate_inserts, -1, true);
        }
        statistics["time- init_inserts"] +=
                duration_cast<duration<double>>(high_resolution_clock::now() -
                                                start_init_inserts)
                        .count();
    }
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    // XGES-0 main loop, in order: reverse, delete, insert; one operator per iteration
    while (!candidate_inserts.empty() || !candidate_reverses.empty() ||
           !candidate_deletes.empty()) {
        edge_modifications.clear();

        if (!candidate_deletes.empty()) {
            // apply the best delete if possible
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();
            if (pdag.is_delete_valid(best_delete)) {
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("{}: {}", i_operations, best_delete);
            } else {
                continue;
            }
        } else if (!candidate_reverses.empty()) {
            // apply the best reverse if possible (no delete available)
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();
            if (pdag.is_reverse_valid(best_reverse, unblocked_paths_map)) {
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                _logger->debug("{}: {}", i_operations, best_reverse);
            } else {
                continue;
            }
        } else if (!candidate_inserts.empty()) {
            // apply the best insert if possible (no delete or reverse available)
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();
            if (best_insert.y == last_insert.y &&
                abs(best_insert.score - last_insert.score) < 1e-10 &&
                best_insert.x == last_insert.x && best_insert.T == last_insert.T) {
                statistics["probable_insert_duplicates"] += 1;
                continue;
            }
            last_insert = std::move(best_insert);
            if (pdag.is_insert_valid(last_insert, unblocked_paths_map)) {
                pdag.apply_insert(last_insert, edge_modifications);
                total_score += last_insert.score;
                _logger->debug("{}: {}", i_operations, last_insert);
            } else {
                continue;
            }
        }
        // here, we have applied an operator
        i_operations++;
        for (auto &edge_modification: edge_modifications) {
            _logger->trace("\tEdge {}", edge_modification.second);
        }
        // update the new possible operators

        // update_operator_candidates_naive(candidate_inserts, candidate_reverses,
        //                                  candidate_deletes);

        // update_operator_candidates(edge_modifications, candidate_inserts,
        //                            candidate_reverses, candidate_deletes);

        update_operator_candidates_v2(edge_modifications, candidate_inserts,
                                      candidate_reverses, candidate_deletes,
                                      unblocked_paths_map);
    }
}

void XGES::update_operator_candidates_naive(std::vector<Insert> &candidate_inserts,
                                            std::vector<Reverse> &candidate_reverses,
                                            std::vector<Delete> &candidate_deletes) {
    candidate_inserts.clear();
    candidate_reverses.clear();
    candidate_deletes.clear();
    for (int y = 0; y < n_variables; ++y) {
        find_inserts_to_y(y, candidate_inserts);
        find_reverse_to_y(y, candidate_reverses);
        find_deletes_to_y(y, candidate_deletes);
    }
}


void XGES::update_operator_candidates(EdgeModificationsMap &edge_modifications,
                                      std::vector<Insert> &candidate_inserts,
                                      std::vector<Reverse> &candidate_reverses,
                                      std::vector<Delete> &candidate_deletes) {
    std::set<int> touched_nodes;
    std::set<int> full_insert_to_y;

    for (auto &edge_modification_key_value: edge_modifications) {
        auto &edge = edge_modification_key_value.second;
        touched_nodes.insert(edge.x);
        touched_nodes.insert(edge.y);
        if (edge.is_reverse() || edge.is_new_undirected()) {
            full_insert_to_y.insert(edge.x);
            full_insert_to_y.insert(edge.y);
            std::ranges::set_intersection(
                    pdag.get_neighbors(edge.x), pdag.get_neighbors(edge.y),
                    std::inserter(full_insert_to_y, full_insert_to_y.begin()));

        } else if (edge.is_new_directed()) {
            if (edge.new_type == EdgeType::DIRECTED_TO_Y) {
                full_insert_to_y.insert(edge.y);
            } else {
                full_insert_to_y.insert(edge.x);
            }
            std::ranges::set_intersection(
                    pdag.get_neighbors(edge.x), pdag.get_neighbors(edge.y),
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
    for (const auto node: full_insert_to_y) {
        find_inserts_to_y(node, candidate_inserts);
    }

    for (auto &edge_modification: edge_modifications) {
        // symmetric in x and y
        auto &edge = edge_modification.second;
        for (auto target: pdag.get_neighbors(edge.y)) {
            if (full_insert_to_y.contains(target)) { continue; }
            find_inserts_to_y(target, candidate_inserts, edge.x);
        }
        for (auto target: pdag.get_neighbors(edge.x)) {
            if (full_insert_to_y.contains(target)) { continue; }
            find_inserts_to_y(target, candidate_inserts, edge.y);
        }

        if (edge.is_new_directed()) {
            int x_ = edge.get_new_source();
            int y_ = edge.get_new_target();
            if (!full_insert_to_y.contains(x_)) {
                for (auto source: pdag.get_adjacent(y_)) {
                    find_inserts_to_y(x_, candidate_inserts, source);
                }
            }
        }
    }
    statistics["time- update_operators"] +=
            duration_cast<duration<double>>(high_resolution_clock::now() - start_time)
                    .count();
}

void add_pairs(std::set<std::pair<int, int>> &pairs, const FlatSet &xs, int y) {
    for (auto x: xs) { pairs.emplace(x, y); }
}
void add_pairs(std::set<std::pair<int, int>> &pairs, int x, const FlatSet &ys) {
    for (auto y: ys) { pairs.emplace(x, y); }
}
void add_pairs(std::set<std::pair<int, int>> &pairs, const FlatSet &xs,
               const FlatSet &ys) {
    for (auto x: xs) {
        for (auto y: ys) { pairs.emplace(x, y); }
    }
}

void XGES::update_operator_candidates_v2(EdgeModificationsMap &edge_modifications,
                                         std::vector<Insert> &candidate_inserts,
                                         std::vector<Reverse> &candidate_reverses,
                                         std::vector<Delete> &candidate_deletes,
                                         UnblockedPathsMap &unblocked_paths_map) {

    // First, undo all the edge modifications
    for (auto &[fst, edge_modification]: edge_modifications) {
        pdag.apply_edge_modification(edge_modification, true);
    }

    std::set<int> full_insert_to_y;
    std::map<int, std::set<int>> partial_insert_to_y;
    std::set<int> full_delete_to_y;
    std::set<int> full_delete_from_x;
    std::set<std::pair<int, int>> delete_x_y;
    std::set<int> full_reverse_to_y;
    std::set<int> full_reverse_from_x;
    std::set<std::pair<int, int>> reverse_x_y;

    // Re-apply the edge modifications and update the operators
    for (auto &[fst, edge_modification]: edge_modifications) {
        int a;
        int b;
        if (edge_modification.is_old_directed()) {
            a = edge_modification.get_old_source();
            b = edge_modification.get_old_target();
        } else if (edge_modification.is_new_directed()) {
            a = edge_modification.get_new_source();
            b = edge_modification.get_new_target();
        } else {
            a = edge_modification.x;
            b = edge_modification.y;
        }
        // Track inserts
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                // y = a
                full_insert_to_y.insert(a);
            case 2:// a  b becomes a → b
                // y = b
                full_insert_to_y.insert(b);
                // y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_insert_to_y, full_insert_to_y.begin()));
                // x=a and y \in Ne(b)
                for (auto target: pdag.get_neighbors(b)) {
                    partial_insert_to_y[target].insert(a);
                }
                // x=b and y \in Ne(a)
                for (auto target: pdag.get_neighbors(a)) {
                    partial_insert_to_y[target].insert(b);
                }
                break;
            case 3:// a -- b becomes a  b
                // x=a and y \in Ne(b) u {b}
                for (auto target: pdag.get_neighbors(b)) {
                    if (target == a) { continue; }
                    partial_insert_to_y[target].insert(a);
                }
                partial_insert_to_y[b].insert(a);
                // y=a and x \in Ad(b)
                partial_insert_to_y[a].insert(pdag.get_adjacent(b).begin(),
                                              pdag.get_adjacent(b).end());
                // x=b and y \in Ne(a) u {a}
                for (auto target: pdag.get_neighbors(a)) {
                    if (target == b) { continue; }
                    partial_insert_to_y[target].insert(b);
                }
                partial_insert_to_y[a].insert(b);
                // y=b and x \in Ad(a)
                partial_insert_to_y[b].insert(pdag.get_adjacent(a).begin(),
                                              pdag.get_adjacent(a).end());
                //SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b});
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a});
                break;
            case 4:// a -- b becomes a → b
                // y = a and x \in Ad(b)
                partial_insert_to_y[a].insert(pdag.get_adjacent(b).begin(),
                                              pdag.get_adjacent(b).end());
                // y = b
                full_insert_to_y.insert(b);
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a});
                break;
            case 5:// a → b becomes a  b
                    // x=a and y \in Ne(b) u {b}
                for (auto target: pdag.get_neighbors(b)) {
                    partial_insert_to_y[target].insert(a);
                }
                partial_insert_to_y[b].insert(a);
                // x=b and y \in Ne(a) u {a}
                for (auto target: pdag.get_neighbors(a)) {
                    partial_insert_to_y[target].insert(b);
                }
                partial_insert_to_y[a].insert(b);
                full_insert_to_y.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b});
                break;
            case 6:// a → b becomes a -- b
                full_insert_to_y.insert(a);
                full_insert_to_y.insert(b);
                break;
            case 7:// a → b becomes a ← b
                full_insert_to_y.insert(a);
                full_insert_to_y.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b});
                break;
            default:
                throw std::runtime_error("Invalid modification");
        }
        // Track deletes
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                full_delete_to_y.insert(a);
                full_delete_to_y.insert(b);
                full_delete_from_x.insert(a);
                full_delete_from_x.insert(b);
                // todo: x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b)

                break;
            case 2:// a  b becomes a → b
                full_delete_to_y.insert(b);
                full_delete_from_x.insert(a);
                full_delete_from_x.insert(b);
                // todo: x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b)
                break;
            case 3:// a -- b becomes a  b
                break;
            case 4:// a -- b becomes a → b
                full_delete_to_y.insert(b);
                break;
            case 5:// a → b becomes a  b
                full_delete_to_y.insert(b);
                break;
            case 6:// a → b becomes a -- b
            case 7:// a → b becomes a ← b
                full_delete_to_y.insert(a);
                full_delete_to_y.insert(b);
                break;
        }

        // Track reverse
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                // y \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                //y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_reverse_to_y, full_reverse_to_y.begin()));
                // x \in {a, b}
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                break;
            case 2:// a  b becomes a → b
                // y = b
                full_reverse_to_y.insert(b);
                // y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_reverse_to_y, full_reverse_to_y.begin()));
                // x \in {a, b}
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                break;
            case 3:// a -- b becomes a  b
                // y \in {a, b} or x \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({b, a});
                }
                break;
            case 4:// a -- b becomes a → b
                // y \in {a, b} or x = b
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(b);
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({b, a});
                }
                break;
            case 5:// a → b becomes a  b
                // y = b or x \in {a, b}
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                break;
            case 6:// a → b becomes a -- b
                // y \in {a, b} or x = b
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(b);
                break;
            case 7:// a → b becomes a ← b
                // y \in {a, b} or x \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                break;
        }
        pdag.apply_edge_modification(edge_modification);
    }

    auto start_time = high_resolution_clock::now();

    // Find the inserts
    // step 1: remove partial inserts that are now full
    std::vector<int> keys_to_erase;
    for (auto [y, xs]: partial_insert_to_y) {
        if (full_insert_to_y.contains(y)) { keys_to_erase.push_back(y); }
    }
    for (auto key: keys_to_erase) { partial_insert_to_y.erase(key); }
    // step 2: find the partial inserts
    for (auto [y, xs]: partial_insert_to_y) {
        for (auto x: xs) {
            // check that x is not adjacent to y
            if (pdag.get_adjacent(y).contains(x)) { continue; }
            find_inserts_to_y(y, candidate_inserts, x);
        }
    }
    // step 3: find the full inserts
    for (auto y: full_insert_to_y) { find_inserts_to_y(y, candidate_inserts); }

    // Find the deletes
    // step 1: form the edges to delete
    for (auto x: full_delete_from_x) { add_pairs(delete_x_y, x, pdag.get_neighbors(x)); }
    for (auto x: full_delete_from_x) { add_pairs(delete_x_y, x, pdag.get_children(x)); }
    for (auto y: full_delete_to_y) { add_pairs(delete_x_y, pdag.get_parents(y), y); }
    for (auto y: full_delete_to_y) { add_pairs(delete_x_y, pdag.get_neighbors(y), y); }
    // step 2: find the deletes
    for (auto [x, y]: delete_x_y) { find_delete_to_y_from_x(y, x, candidate_deletes); }

    // Find the reverses
    // step 1: form the edges to reverse
    for (auto x: full_reverse_from_x) { add_pairs(reverse_x_y, x, pdag.get_parents(x)); }
    for (auto y: full_reverse_to_y) { add_pairs(reverse_x_y, pdag.get_children(y), y); }
    // step 2: find the reverses
    for (auto [x, y]: reverse_x_y) {
        // check that x ← y
        if (!pdag.has_directed_edge(y, x)) { continue; }
        find_reverse_to_y_from_x(y, x, candidate_reverses);
    }


    statistics["time- update_operators"] +=
            duration_cast<duration<double>>(high_resolution_clock::now() - start_time)
                    .count();
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
 * @param parent_x
 * @param positive_only
 */

void XGES::find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts, int parent_x,
                             bool positive_only) {
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y = pdag.get_parents(y);

    std::set<int> possible_parents;

    if (parent_x != -1) {
        possible_parents.insert(parent_x);
    } else {
        // for now: no pre-selection
        auto &nodes = pdag.get_nodes_variables();

        // 1. x is not adjacent to y (x ∉ Ad(y))
        std::ranges::set_difference(
                nodes, adjacent_y,
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
        FlatSet &effective_parents_y =
                neighbors_y_adjacent_x;// just renaming it, no copy necessary
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
                    duration_cast<duration<double>>(high_resolution_clock::now() - start)
                            .count();
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
                if (std::ranges::includes(adjacent_z, T) &&
                    std::ranges::includes(adjacent_z, neighbors_y_adjacent_x)) {
                    // T' is a candidate
                    auto T_prime = T;
                    T_prime.insert(z);
                    auto effective_parents_prime = effective_parents;
                    effective_parents_prime.insert(z);
                    stack.emplace(std::move(T_prime), it,
                                  std::move(effective_parents_prime));
                }
            }
        }
    }
}

void XGES::find_delete_to_y_from_x(int y, int x, std::vector<Delete> &candidate_deletes,
                                   bool positive_only) const {
    const FlatSet &parents_y = pdag.get_parents(y);
    auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
    bool directed_xy = pdag.has_directed_edge(x, y);

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
        if (score > 0 || !positive_only) {
            candidate_deletes.emplace_back(x, y, O, score, effective_parents,
                                           directed_xy);
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
            if (std::ranges::includes(adjacent_z, O)) {
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
 * @param positive_only
 */
void XGES::find_deletes_to_y(const int y, std::vector<Delete> &candidate_deletes,
                             bool positive_only) const {
    auto &neighbors_y = pdag.get_neighbors(y);
    auto &parents_y = pdag.get_parents(y);

    for (int x: parents_y) {
        find_delete_to_y_from_x(y, x, candidate_deletes, positive_only);
    }
    for (int x: neighbors_y) {
        find_delete_to_y_from_x(y, x, candidate_deletes, positive_only);
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
 * @param candidate_reverses
 */

// can i just do insert x -> y?
void XGES::find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses) {
    // look for all possible x ← y
    auto &children_y = pdag.get_children(y);

    for (int x: children_y) {
        auto &parents_x = pdag.get_parents(x);
        std::vector<Insert> candidate_inserts;
        find_inserts_to_y(y, candidate_inserts, x, false);

        for (auto &insert: candidate_inserts) {
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
        find_inserts_to_y(y, candidate_inserts, x, false);

        for (auto &insert: candidate_inserts) {
            // change if we parallelize
            double score = insert.score + scorer->score_delete(x, parents_x, y);

            if (score > 0) {
                candidate_reverses.emplace_back(insert, score, parents_x);
                std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
            }
        }
    }
}

void XGES::find_reverse_to_y_from_x(int y, int x,
                                    std::vector<Reverse> &candidate_reverses) {
    if (!pdag.has_directed_edge(y, x)) { return; }
    std::vector<Insert> candidate_inserts;
    find_inserts_to_y(y, candidate_inserts, x, false);
    auto &parents_x = pdag.get_parents(x);
    for (auto &insert: candidate_inserts) {
        double score = insert.score + scorer->score_delete(x, parents_x, y);
        if (score > 0) {
            candidate_reverses.emplace_back(insert, score, parents_x);
            std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
        }
    }
}


const PDAG &XGES::get_pdag() const { return pdag; }

double XGES::get_score() const { return total_score; }
double XGES::get_initial_score() const{ return initial_score; }
