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

/** @brief Constructor of the XGES class.
 *
 * @param n_variables The number of variables in the dataset.
 * @param scorer The scorer to use to score the PDAGs.
 *
 */
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


/** @brief Fit the XGES algorithm using the scorer provided in the constructor.
 *
 * @param extended_search If true, the extended search of XGES is performed. Otherwise,
 * only XGES-0 is performed.
 */
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

/** @brief Extended search of XGES.
 *
 * Extended search is performed after XGES-0 has exhausted all possible operators.
 * It successively delete some edges, and resume the search with XGES-0 (preventing
 * the deleted edges from being inserted again). If the score increases, the new PDAG
 * is kept.
 *
 * At the beggining of the search, `all_edge_deletes` contains all the possible deletes
 * for the current PDAG. If the PDAG is updated, we keep iterating through `all_edge_deletes`
 * until all the deletes have been tried. Only then we recompute `all_edge_deletes` on the
 * last updated PDAG. Until no new deletes are found. This is more efficient than recomputing
 * all the deletes at each PDAG update, and re-testing all of them (as most of them will be
 * the same).
 *
 *
 * @param unblocked_paths_map The unblocked paths map used by XGES-0.
 */
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

        // xges_copy.update_operator_candidates_naive(candidate_inserts,
        // candidate_reverses, candidate_deletes);
        xges_copy.update_operator_candidates_efficient(
                edge_modifications, candidate_inserts, candidate_reverses,
                candidate_deletes, blocked_paths_map_copy);
        xges_copy.heuristic_xges0(candidate_inserts, candidate_reverses,
                                  candidate_deletes, blocked_paths_map_copy, false);
        if (pdag == xges_copy.pdag) { continue; }
        if (xges_copy.total_score - total_score > 1e-7) {
            _logger->debug("EXTENDED SEARCH ACCEPTED: with increase {} and {}",
                           xges_copy.total_score - total_score, delete_);
            total_score = xges_copy.total_score;
            pdag = xges_copy.pdag;
            unblocked_paths_map = std::move(blocked_paths_map_copy);
            deletes_of_pdag_are_updated = false;
            statistics["extended_search-accepted"] += 1;
        } else {
            _logger->debug("EXTENDED SEARCH REJECTED: {} {}", delete_,
                           xges_copy.total_score);
            statistics["extended_search-rejected"] += 1;
        }
    }
}

/** @brief Run the XGES-0 heuristic.
 *
 * XGES-0 is the main heuristic of the XGES algorithm. It maintains lists of all possible
 * operators and apply them in order: delete, reverse, and insert. XGES-0 stops when no
 * more operators can be applied.
 * XGES-0 uses update_operator_candidates_efficient to update the operators after each
 * PDAG update.
 *
 * @param candidate_inserts The vector of candidate inserts.
 * @param candidate_reverses The vector of candidate reverses.
 * @param candidate_deletes The vector of candidate deletes.
 * @param unblocked_paths_map The unblocked paths map used by XGES-0.
 * @param initialize_inserts If true, the candidate inserts are initialized at the beginning
 * of the heuristic.
 */
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
        statistics["initialize_inserts-time"] += measure_time(start_init_inserts);
    }
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    // XGES-0 main loop, in order: delete, reverse, insert; one operator per iteration
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
        // if we reach this point, we have applied an operator
        i_operations++;
        for (auto &edge_modification: edge_modifications) {
            _logger->trace("\tEdge {}", edge_modification.second);
        }

        // update the new possible operators
        // update_operator_candidates_naive(candidate_inserts, candidate_reverses,
        //                                  candidate_deletes);
        update_operator_candidates_efficient(edge_modifications, candidate_inserts,
                                             candidate_reverses, candidate_deletes,
                                             unblocked_paths_map);
    }
}

void XGES::fit_ops(bool use_reverse) {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_inserts.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    for (int y = 0; y < n_variables; ++y) {
        find_inserts_to_y(y, candidate_inserts, -1, true);
    }
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    while (!candidate_inserts.empty() || !candidate_reverses.empty() ||
           !candidate_deletes.empty()) {

        // Look for the operator with the highest score
        // Each of candidate_inserts, candidate_reverses, candidate_deletes are heaps
        // with the highest score at the front.
        double best_score = 0;
        int best_type = -1;// 0: insert, 1: reverse, 2: delete

        if (!candidate_inserts.empty()) {
            if (candidate_inserts.front().score > best_score) {
                best_type = 0;
                best_score = candidate_inserts.front().score;
            }
        }
        if (!candidate_reverses.empty() && use_reverse) {
            if (candidate_reverses.front().score > best_score) {
                best_type = 1;
                best_score = candidate_reverses.front().score;
            }
        }
        if (!candidate_deletes.empty()) {
            if (candidate_deletes.front().score > best_score) {
                best_type = 2;
                best_score = candidate_deletes.front().score;
            }
        }

        edge_modifications.clear();
        if (best_type == 0) {
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
        } else if (best_type == 1) {
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
        } else if (best_type == 2) {
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
        } else {
            continue;
        }


        i_operations++;
        for (auto &edge_modification: edge_modifications) {
            _logger->trace("\tEdge {}", edge_modification.second);
        }

        update_operator_candidates_naive(candidate_inserts, candidate_reverses,
                                         candidate_deletes);

        if (!use_reverse) { candidate_reverses.clear(); }
    }
}

void XGES::fit_ges(bool use_reverse) {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_inserts.reserve(n_variables);
    candidate_deletes.reserve(n_variables);


    EdgeModificationsMap edge_modifications;
    int i_operations = 1;
    int last_operation = -1;

    while (i_operations > last_operation) {
        last_operation = i_operations;

        // only apply inserts
        for (int y = 0; y < n_variables; ++y) { find_inserts_to_y(y, candidate_inserts); }
        while (!candidate_inserts.empty()) {
            edge_modifications.clear();
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();
            if (pdag.is_insert_valid(best_insert, unblocked_paths_map)) {
                pdag.apply_insert(best_insert, edge_modifications);
                total_score += best_insert.score;
                _logger->debug("{}: {}", i_operations, best_insert);
                i_operations++;
                candidate_inserts.clear();
                for (int y = 0; y < n_variables; ++y) {
                    find_inserts_to_y(y, candidate_inserts);
                }
            }
        }
        // only apply deletes
        for (int y = 0; y < n_variables; ++y) { find_deletes_to_y(y, candidate_deletes); }
        while (!candidate_deletes.empty()) {
            edge_modifications.clear();
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();
            if (pdag.is_delete_valid(best_delete)) {
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("{}: {}", i_operations, best_delete);
                i_operations++;
                candidate_deletes.clear();
                for (int y = 0; y < n_variables; ++y) {
                    find_deletes_to_y(y, candidate_deletes);
                }
            }
        }
        // only apply reverses, if use_reverse is true
        // note that without reverse, there is no multiple
        // iterations of forward and backward
        if (!use_reverse) { return; }
        for (int y = 0; y < n_variables; ++y) {
            find_reverse_to_y(y, candidate_reverses);
        }
        while (!candidate_reverses.empty()) {
            edge_modifications.clear();
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();
            if (pdag.is_reverse_valid(best_reverse, unblocked_paths_map)) {
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                _logger->debug("{}: {}", i_operations, best_reverse);
                i_operations++;
                candidate_reverses.clear();
                for (int y = 0; y < n_variables; ++y) {
                    find_reverse_to_y(y, candidate_reverses);
                }
            }
        }
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

/** @brief Update the candidate operators after a PDAG update.
 *
 * After a PDAG update with a set of edge modifications, the candidate operators
 * are updated following the conditions detailed in the XGES paper. The insert, reverse,
 * and delete operators are updated according to the type of each edge modification.
 *
 * @param edge_modifications The edge modifications that have been applied to the PDAG.
 * @param candidate_inserts The vector of candidate inserts.
 * @param candidate_reverses The vector of candidate reverses.
 * @param candidate_deletes The vector of candidate deletes.
 * @param unblocked_paths_map The unblocked paths map used by XGES-0.
 */
void XGES::update_operator_candidates_efficient(EdgeModificationsMap &edge_modifications,
                                                std::vector<Insert> &candidate_inserts,
                                                std::vector<Reverse> &candidate_reverses,
                                                std::vector<Delete> &candidate_deletes,
                                                UnblockedPathsMap &unblocked_paths_map) {

    auto start_time = high_resolution_clock::now();
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

    // Re-apply the edge modifications one by one and update the operators
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
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a}); // Leave them for the reverse
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
                // unblocked_paths_map.erase({b, a}); // Leave them for the reverse
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
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
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
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                break;
        }
        // Track deletes
        FlatSet x_intersection;
        FlatSet y_intersection;
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                full_delete_to_y.insert(a);
                full_delete_to_y.insert(b);
                full_delete_from_x.insert(a);
                full_delete_from_x.insert(b);

                // x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b) [almost never happens]
                std::ranges::set_intersection(
                        pdag.get_adjacent(a), pdag.get_adjacent(b),
                        std::inserter(x_intersection, x_intersection.begin()));
                if (!x_intersection.empty()) {
                    std::ranges::set_intersection(
                            pdag.get_neighbors(a), pdag.get_neighbors(b),
                            std::inserter(y_intersection, y_intersection.begin()));
                    add_pairs(delete_x_y, x_intersection, y_intersection);
                }

                break;
            case 2:// a  b becomes a → b
                // TODO: merge (appropriately) with case 1
                full_delete_to_y.insert(b);
                full_delete_from_x.insert(a);
                full_delete_from_x.insert(b);
                // x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b) [almost never happens]
                std::ranges::set_intersection(
                        pdag.get_adjacent(a), pdag.get_adjacent(b),
                        std::inserter(x_intersection, x_intersection.begin()));
                if (!x_intersection.empty()) {
                    std::ranges::set_intersection(
                            pdag.get_neighbors(a), pdag.get_neighbors(b),
                            std::inserter(y_intersection, y_intersection.begin()));
                    add_pairs(delete_x_y, x_intersection, y_intersection);
                }
                break;
            case 3:// a -- b becomes a  b
                break;
            case 4:// a -- b becomes a → b
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
            if (!pdag.get_adjacent(y).contains(x) && x != y) {
                find_inserts_to_y(y, candidate_inserts, x);
            }
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
    for (auto [x, y]: delete_x_y) {
        if (x != y) { find_delete_to_y_from_x(y, x, candidate_deletes); }
    }

    // Find the reverses
    // step 1: form the edges to reverse
    for (auto x: full_reverse_from_x) { add_pairs(reverse_x_y, x, pdag.get_parents(x)); }
    for (auto y: full_reverse_to_y) { add_pairs(reverse_x_y, pdag.get_children(y), y); }
    // step 2: find the reverses
    for (auto [x, y]: reverse_x_y) {
        // check that x ← y
        if (pdag.has_directed_edge(y, x) && x != y) {
            find_reverse_to_y_from_x(y, x, candidate_reverses);
        }
    }
    statistics["update_operators-time"] += measure_time(start_time);
}


/**
 * @brief Find and score all possible inserts to y.
 *
 * The candidate  Insert(x, y, T, E) are such that:
 *  1. x is not adjacent to y (x ∉ Ad(y))
 *  2. T ⊆ Ne(y) \ Ad(x)
 *  3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 *  Not enforced at that stage: [Ne(y) ∩ Ad(x)] ∪ T blocks all semi-directed paths from y to x
 *  5. E = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
 *
 * @param y The target node.
 * @param candidate_inserts The vector of candidate inserts to append inserts to.
 * @param parent_x The parent of y, if known. If -1, no pre-selection is done.
 * @param positive_only If true, only inserts with a positive score are considered.
 */
void XGES::find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts, int parent_x,
                             bool positive_only) const {
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
            double score = scorer->score_insert(y, effective_parents, x);
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

/** Find all possible deletes to y from x.
 *
 * The candidate deletes (x, y, C, E) are such that:
 *  1. x is a parent or a neighbor of y
 *  2. C ⊆ [Ne(y) ∩ Ad(x)]
 *  3. C is a clique
 *  4. E = Pa(y) ∪ [Ne(y) ∩ Ad(x)] \ C ∪ {x}
 *
 *
 * @param y The target node.
 * @param x The source node.
 * @param candidate_deletes The vector of candidate deletes to append deletes to.
 * @param positive_only If true, only deletes with a positive score are considered.
 */
void XGES::find_delete_to_y_from_x(int y, int x, std::vector<Delete> &candidate_deletes,
                                   bool positive_only) const {
    const FlatSet &parents_y = pdag.get_parents(y);
    auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
    bool directed_xy = pdag.has_directed_edge(x, y);

    // find all possible C ⊆ [Ne(y) ∩ Ad(x)] that are cliques
    // <C set of nodes, iterator over neighbors_y_adjacent_x, set of effective_parents>
    std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
    FlatSet effective_parents_init;
    effective_parents_init.reserve(parents_y.size() + neighbors_y_adjacent_x.size() + 1);
    // note: Chickering Corollary 18 is incorrect. Pa(y) might not contain x, it has to be added.
    union_with_single_element(parents_y, x, effective_parents_init);
    // note that C = {} is valid
    stack.emplace(FlatSet{}, neighbors_y_adjacent_x.begin(), effective_parents_init);

    while (!stack.empty()) {
        auto top = std::move(stack.top());
        stack.pop();
        auto C = std::get<0>(top);
        auto it = std::get<1>(top);
        auto effective_parents = std::get<2>(top);

        double score = scorer->score_delete(y, effective_parents, x);
        if (score > 0 || !positive_only) {
            candidate_deletes.emplace_back(x, y, C, score, effective_parents,
                                           directed_xy);
            std::push_heap(candidate_deletes.begin(), candidate_deletes.end());
        }

        // Look for other candidate C using the iterator, which gives us the next elements
        // z to consider.
        while (it != neighbors_y_adjacent_x.end()) {
            // We define C' = C ∪ {z} and we check if C' is a clique.
            // Since C is a clique, we only need to check that z is adjacent to all nodes
            // in C.
            auto z = *it;
            ++it;
            auto &adjacent_z = pdag.get_adjacent(z);
            // We check C ⊆ Ad(z)
            if (std::ranges::includes(adjacent_z, C)) {
                // C' is a candidate
                auto C_prime = C;
                C_prime.insert(z);
                auto effective_parents_prime = effective_parents;
                effective_parents_prime.insert(z);
                stack.emplace(std::move(C_prime), it, std::move(effective_parents_prime));
            }
        }
    }
}

/** @brief Find all possible deletes to y.
 *
 * @see find_delete_to_y_from_x
 *
 * @param y The target node.
 * @param candidate_deletes The vector of candidate deletes to append deletes to.
 * @param positive_only If true, only deletes with a positive score are considered.
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
 * Only used for the naive update of the operators.
 */
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

/** @brief Find all possible reverses to y from x.
 *
 * @see find_inserts_to_y
 *
 * The candidate Reverse(x, y, T, E, F) are such that:
 * 1. x ← y
 * 2. T ⊆ Ne(y) \ Ad(x)
 * 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 * 4. All semi-directed paths from y to x are blocked by [Ne(y) ∩ Ad(x)] ∪ T ∪ Ne(x)
 * 5. E = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
 * 6. F = Pa(x)
 *
 * A reverse is composed of a delete and an insert. The function first finds all possible
 * inserts to y from x using `find_inserts_to_y` and then evaluate the delete part.
 *
 */
void XGES::find_reverse_to_y_from_x(int y, int x,
                                    std::vector<Reverse> &candidate_reverses) const {
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
