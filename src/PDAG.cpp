//
// Created by Achille Nazaret on 11/3/23.
//
#include "PDAG.h"
#include "set_ops.h"
#include "utils.h"
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace std::chrono;

PDAG::PDAG(const int num_nodes)
    : num_variables(num_nodes), _block_semi_directed_path_queue(num_nodes) {
    for (int i = 0; i < num_nodes; i++) {
        nodes_variables.push_back(i);
        nodes_all.push_back(i);
        children.emplace_back();
        parents.emplace_back();
        neighbors.emplace_back();
        adjacent.emplace_back();
        adjacent_reachable.emplace_back();

        forbidden_insert_parents.emplace_back();
        forbidden_insert_children.emplace_back();
    }
    _block_semi_directed_path_visited.resize(num_nodes);
    _block_semi_directed_path_blocked.resize(num_nodes);
    _block_semi_directed_path_parent.resize(num_nodes);
}

int PDAG::get_number_of_edges() const {
    return number_of_directed_edges + number_of_undirected_edges;
}

const std::vector<int> &PDAG::get_nodes_variables() const { return nodes_variables; }

std::vector<std::pair<int, int>> PDAG::get_directed_edges() const {
    std::vector<std::pair<int, int>> edges;
    for (const auto &node: nodes_variables) {
        for (int child: children.at(node)) { edges.emplace_back(node, child); }
    }
    return edges;
}

std::vector<std::pair<int, int>> PDAG::get_undirected_edges() const {
    std::vector<std::pair<int, int>> edges;
    for (const auto &node: nodes_variables) {
        for (int neighbor: neighbors.at(node)) {
            if (node < neighbor) { edges.emplace_back(node, neighbor); }
        }
    }
    return edges;
}

const FlatSet &PDAG::get_parents(const int node) const { return parents[node]; }

const FlatSet &PDAG::get_children(const int node) const { return children[node]; }

const FlatSet &PDAG::get_neighbors(const int node) const { return neighbors[node]; }

const FlatSet &PDAG::get_adjacent(const int node) const { return adjacent[node]; }

const FlatSet &PDAG::get_adjacent_reachable(const int node) const {
    return adjacent_reachable.at(node);
}

FlatSet PDAG::get_neighbors_adjacent(const int node_y, const int node_x) const {
    FlatSet result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);
    std::ranges::set_intersection(neighbors_set, adjacent_set,
                                  std::inserter(result, result.begin()));
    return result;
}

FlatSet PDAG::get_neighbors_not_adjacent(const int node_y, const int node_x) const {
    FlatSet result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);
    std::ranges::set_difference(neighbors_set, adjacent_set,
                                std::inserter(result, result.begin()));
    return result;
}

bool PDAG::has_directed_edge(int x, int y) const { return children.at(x).contains(y); }

bool PDAG::has_undirected_edge(int x, int y) const { return neighbors.at(x).contains(y); }

void PDAG::remove_directed_edge(int x, int y) {
    children[x].erase(y);
    parents[y].erase(x);
    adjacent[x].erase(y);
    adjacent[y].erase(x);
    adjacent_reachable[x].erase(y);
    number_of_directed_edges--;
}

void PDAG::remove_undirected_edge(int x, int y) {
    neighbors[x].erase(y);
    neighbors[y].erase(x);
    adjacent[x].erase(y);
    adjacent[y].erase(x);
    adjacent_reachable[x].erase(y);
    adjacent_reachable[y].erase(x);
    number_of_undirected_edges--;
}

void PDAG::add_directed_edge(int x, int y) {
    children[x].insert(y);
    parents[y].insert(x);
    adjacent[x].insert(y);
    adjacent[y].insert(x);
    adjacent_reachable[x].insert(y);
    number_of_directed_edges++;
}

void PDAG::add_undirected_edge(int x, int y) {
    neighbors[x].insert(y);
    neighbors[y].insert(x);
    adjacent[x].insert(y);
    adjacent[y].insert(x);
    adjacent_reachable[x].insert(y);
    adjacent_reachable[y].insert(x);
    number_of_undirected_edges++;
}

void PDAG::apply_edge_modification(EdgeModification edge_modification, bool undo) {
    EdgeType old_type = edge_modification.old_type;
    EdgeType new_type = edge_modification.new_type;
    if (undo) { std::swap(old_type, new_type); }
    switch (old_type) {
        case EdgeType::UNDIRECTED:
            remove_undirected_edge(edge_modification.x, edge_modification.y);
            break;
        case EdgeType::DIRECTED_TO_X:
            remove_directed_edge(edge_modification.y, edge_modification.x);
            break;
        case EdgeType::DIRECTED_TO_Y:
            remove_directed_edge(edge_modification.x, edge_modification.y);
            break;
        default:
            break;
    }
    switch (new_type) {
        case EdgeType::UNDIRECTED:
            add_undirected_edge(edge_modification.x, edge_modification.y);
            break;
        case EdgeType::DIRECTED_TO_X:
            add_directed_edge(edge_modification.y, edge_modification.x);
            break;
        case EdgeType::DIRECTED_TO_Y:
            add_directed_edge(edge_modification.x, edge_modification.y);
            break;
        default:
            break;
    }
}


/**
 * @brief Check if the nodes in the set form a clique in the PDAG (i.e., are all adjacent to each other).
 *
 * @param nodes_subset The set of nodes to check if they form a clique
 * @return `true` if the nodes form a clique, `false` otherwise.
 */
bool PDAG::is_clique(const FlatSet &nodes_subset) const {
    for (int node1: nodes_subset) {
        const FlatSet &adjacent_set = get_adjacent(node1);
        for (int node2: nodes_subset) {
            if (node1 != node2 && adjacent_set.find(node2) == adjacent_set.end()) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Check if `blocked_nodes` block all semi-directed paths from `y` to `x` in the PDAG.
 *
 * If `ignore_direct_edge` is `true`, the possible direct edge from `y` to `x` is ignored,
 * which is useful for checking reverse operators.
 *
 * If a non-blocked path is found, each edge e in the path is added to `unblocked_paths_map`,
 * and maps to to the pair (x, y). If in the future, e gets deleted (or oriented right-to-left),
 * the path from y to x will be blocked, and XGES should re-test inserts and reverses from
 * y to x.
 *
 * @param y The source node of the semi-directed paths
 * @param x The destination node of the semi-directed paths
 * @param blocked_nodes The set of nodes that are blocking the paths
 * @param unblocked_paths_map The map of unblocked paths edges
 * @param ignore_direct_edge If the direct edge from y to x should be ignored
 * @return `true` if blocked_nodes block all semi-directed paths from src to dst, `false` otherwise.
 */
bool PDAG::is_blocking_semi_directed_paths(const int y, const int x,
                                           const FlatSet &blocked_nodes,
                                           UnblockedPathsMap &unblocked_paths_map,
                                           const bool ignore_direct_edge) {
    if (y == x) { return false; }
    statistics["block_semi_directed_paths-#calls"] += 1;
    auto start_time = high_resolution_clock::now();

    // Instant lookup for visited and blocked nodes
    auto &visited = _block_semi_directed_path_visited;
    std::ranges::fill(visited, 0);
    auto &blocked = _block_semi_directed_path_blocked;
    std::ranges::fill(blocked, 0);
    for (int n: blocked_nodes) { blocked[n] = 1; }

    // BFS search from y to x, using adjacent_reachable edges, avoiding blocked nodes
    visited[y] = 1;

    auto &queue = _block_semi_directed_path_queue;
    queue.clear();
    queue.push_back(y);

    while (!queue.empty()) {
        int node = queue.pop_front();
        auto &reachable = get_adjacent_reachable(node);

        for (int n: reachable) {
            if (visited[n] || blocked[n] || (node == y && n == x && ignore_direct_edge)) {
                continue;
            }
            _block_semi_directed_path_parent[n] = node;
            if (n == x) {
                statistics["block_semi_directed_paths-false-#"] += 1;
                statistics["block_semi_directed_paths-false-time"] +=
                        measure_time(start_time);
                // retrieve the path
                int current = x;
                while (current != y) {
                    int parent = _block_semi_directed_path_parent[current];
                    unblocked_paths_map[{parent, current}].emplace(x, y);
                    current = parent;
                }
                return false;
            }
            queue.push_back(n);
            visited[n] = 1;
        }
    }
    statistics["block_semi_directed_paths-true-#"] += 1;
    statistics["block_semi_directed_paths-true-time"] += measure_time(start_time);
    return true;
}


/**
 * @brief  Verify if the insert is valid (or that the insert constituting a reverse is valid).
 *
 * Insert(x, y, T, E) is valid if and only if: [in approximate order of complexity, and
 * numbered as in the paper]
 *  0. the edge is not forbidden (the extended search forbid some inserts)
 *  1. x and y are not adjacent (or x is a child of y for reverse)
 *  2. T is a subset of Ne(y) \ Ad(x)
 *  5. E = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
 *  3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 *  4. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x (or
 *  [Ne(y) ∩ Ad(x)] ∪ T ∪ Ne(x) block all semi-directed paths from x to y for reverse)
 *
 * @param insert The insert to verify
 * @param unblocked_paths_map The map of unblocked paths (see is_blocking_semi_directed_paths)
 * @param reverse If the insert is part of a reverse
 * @return `true` if the insert is valid, `false` otherwise.
 */
bool PDAG::is_insert_valid(const Insert &insert, UnblockedPathsMap &unblocked_paths_map,
                           bool reverse) {
    auto start_time = high_resolution_clock::now();
    statistics["is_insert_valid-#calls"] += 1;
    int x = insert.x;
    int y = insert.y;
    auto &T = insert.T;

    // 0. check if edge is forbidden
    if (is_insert_forbidden(x, y)) {
        statistics["is_insert_valid-false_0-#"] += 1;
        statistics["is_insert_valid-false-time"] += measure_time(start_time);
        return false;
    }

    auto &adjacent_x = adjacent.at(x);
    if (!reverse) {
        // 1. x and y are not adjacent
        if (adjacent_x.contains(y)) {
            statistics["is_insert_valid-false_1a-#"] += 1;
            statistics["is_insert_valid-false-time"] += measure_time(start_time);
            return false;
        }
    } else {
        // 1. x ← y
        auto &parents_x = parents.at(x);
        if (!parents_x.contains(y)) {
            statistics["is_insert_valid-false_1b-#"] += 1;
            statistics["is_insert_valid-false-time"] += measure_time(start_time);
            return false;
        }
    }

    // 2. T ⊆ Ne(y) \ Ad(x)
    // <=> T ⊆ Ne(y) and T does not intersect Ad(x)
    auto &neighbors_y = neighbors.at(y);
    if (!is_subset(T, neighbors_y) || have_overlap(T, adjacent_x)) {
        statistics["is_insert_valid-false_2-#"] += 1;
        statistics["is_insert_valid-false-time"] += measure_time(start_time);
        return false;
    }

    // 5. E (insert.effective_parents) == [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    FlatSet ne_y_ad_x_T;
    std::ranges::set_intersection(neighbors_y, adjacent_x,
                                  std::inserter(ne_y_ad_x_T, ne_y_ad_x_T.begin()));
    ne_y_ad_x_T.insert(T.begin(), T.end());
    if (!a_equals_union_b1_b2(insert.effective_parents, ne_y_ad_x_T, parents.at(y))) {
        statistics["is_insert_valid-false_5-#"] += 1;
        statistics["is_insert_valid-false-time"] += measure_time(start_time);
        return false;
    }

    // 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
    if (!is_clique(ne_y_ad_x_T)) {
        statistics["is_insert_valid-false_3-#"] += 1;
        statistics["is_insert_valid-false-time"] += measure_time(start_time);
        return false;
    }

    // 4. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x
    bool ignore_direct_edge = reverse;
    if (reverse) {
        // for reverse: ne_y_ad_x_T is actually [Ne(y) ∩ Ad(x)] ∪ T ∪ Ne(x)
        ne_y_ad_x_T.insert(neighbors.at(x).begin(), neighbors.at(x).end());
    };
    if (!is_blocking_semi_directed_paths(y, x, ne_y_ad_x_T, unblocked_paths_map,
                                         ignore_direct_edge)) {
        statistics["is_insert_valid-false_4-#"] += 1;
        statistics["is_insert_valid-false-time"] += measure_time(start_time);
        return false;
    }

    statistics["is_insert_valid-true-#"] += 1;
    statistics["is_insert_valid-true-time"] += measure_time(start_time);
    return true;
}

/**
 * @brief Verify if the reverse is valid.
 *
 * Reverse(x, y, T, E, F) is valid if and only if:
 * 1. F (parents_x) = Pa(x)
 * 2. the induced insert Insert(x, y, T, E) is valid
 *
 * @param reverse The reverse to verify
 * @param unblocked_paths_map The map of unblocked paths (see is_blocking_semi_directed_paths)
 * @return `true` if the reverse is valid, `false` otherwise.
 */
bool PDAG::is_reverse_valid(const Reverse &reverse,
                            UnblockedPathsMap &unblocked_paths_map) {
    // is Pa(x) unchanged
    if (get_parents(reverse.insert.x) != reverse.parents_x) { return false; }
    return is_insert_valid(reverse.insert, unblocked_paths_map, true);
}

/**
 * @brief Verify if the delete is valid.
 *
 * Delete(x, y, C) is valid if and only if [in approximate order of complexity, and
 * numbered as in the paper]:
 * 1. x and y are neighbors or x is a parent of y
 * 2. C is a subset of Ne(y) ∩ Ad(x)
 * 4. E = C ∪ Pa(y) ∪ {x}
 * 3. C is a clique
 *
 * @param delet The delete to verify
 * @return `true` if the delete is valid, `false` otherwise.
 */
bool PDAG::is_delete_valid(const Delete &delet) const {
    // 1. x and y are neighbors or x is a parent of y [aka y is adjacent_reachable from x]
    const int x = delet.x;
    const int y = delet.y;
    if (!adjacent_reachable.at(x).contains(y)) { return false; }

    // 2. C is a subset of Ne(y) ∩ Ad(x)
    // <=> C ⊆ Ne(y) and C ⊆ Ad(x)
    auto &neighbors_y = neighbors.at(y);
    auto &adjacent_x = adjacent.at(x);
    if (!is_subset(delet.C, neighbors_y) || !is_subset(delet.C, adjacent_x)) {
        return false;
    }

    // 3. E (delet.effective_parents) = C ∪ Pa(y) ∪ {x}
    if (!a_equals_union_b1_b2_and_singleton(delet.effective_parents, delet.C,
                                            parents.at(y), x)) {
        return false;
    }

    // 4. C is a clique
    if (!is_clique(delet.C)) { return false; }

    return true;
}

/**
 * @brief Apply the insert to the PDAG, and complete the resulting PDAG into a CPDAG.
 *
 * Insert(x, y, T, E) does:
 *  1. insert the directed edge x → y
 *  2. for each t ∈ T: orient the (previously undirected) edge between t and y as t → y
 *
 *  The modified edges are added to the edge_modifications_map.
 *
 * @param insert The insert to apply
 * @param edge_modifications_map The map of edge modifications
 */
void PDAG::apply_insert(const Insert &insert,
                        EdgeModificationsMap &edge_modifications_map) {
    const auto start_time = high_resolution_clock::now();
    int x = insert.x;
    int y = insert.y;
    auto &T = insert.T;

    EdgeQueueSet edges_to_check;

    // a. insert the directed edge x → y
    add_directed_edge(x, y);
    edge_modifications_map.update_edge_directed(x, y, EdgeType::NONE);
    // b. for each t ∈ T: orient the (previously undirected) edge between t and y as t → y
    for (int t: T) {
        remove_undirected_edge(t, y);
        add_directed_edge(t, y);
        edge_modifications_map.update_edge_directed(t, y, EdgeType::UNDIRECTED);
    }

    // check if the orientation of the surrounding edges should be updated
    edges_to_check.push_directed(x, y);
    add_adjacent_edges(x, y, edges_to_check);
    // edges t → y are part of a v-structure with x, so we don't need to check them
    for (int t: T) { add_adjacent_edges(t, y, edges_to_check); }

    complete_cpdag(edges_to_check, edge_modifications_map);
    statistics["apply_insert-time"] += measure_time(start_time);
}

/**
 * @brief Apply the reverse to the PDAG, and complete the resulting PDAG into a CPDAG.
 *
 * Reverse(x, y, T, E, F) does:
 *  1. remove the directed edge x ← y
 *  2. apply the insert Insert(x, y, T, E)
 *
 *  The modified edges are added to the edge_modifications_map.
 *
 * @param reverse The reverse to apply
 * @param edge_modifications_map The map of edge modifications
 */
void PDAG::apply_reverse(const Reverse &reverse,
                         EdgeModificationsMap &edge_modifications_map) {
    int x = reverse.insert.x;
    int y = reverse.insert.y;
    // 1. remove the directed edge y → x
    remove_directed_edge(y, x);
    edge_modifications_map.update_edge_none(x, y, EdgeType::DIRECTED_TO_X);
    // 2. apply the insert
    apply_insert(reverse.insert, edge_modifications_map);
}

/**
 * @brief Apply the delete to the PDAG, and complete the resulting PDAG into a CPDAG.
 *
 * Delete(x, y, C) does:
 *  1. remove the directed edge x → y or the undirected edge x - y
 *  2. compute H = Ne(y) ∩ Ad(x) \ C, then for each h ∈ H:
 *    - orient the previously undirected edges between h and y as y → h
 *    - orient any previously undirected edges between x and h as x → h
 *
 *  The modified edges are added to the edge_modifications_map.
 *
 * @param delet The delete to apply
 * @param edge_modifications_map The map of edge modifications
 */
void PDAG::apply_delete(const Delete &delet,
                        EdgeModificationsMap &edge_modifications_map) {
    auto start_time = high_resolution_clock::now();
    if (has_directed_edge(delet.x, delet.y)) {
        // 1. remove the directed edge x → y
        remove_directed_edge(delet.x, delet.y);
        edge_modifications_map.update_edge_none(delet.x, delet.y,
                                                EdgeType::DIRECTED_TO_Y);
    } else {
        // 1. remove the undirected edge x - y
        remove_undirected_edge(delet.x, delet.y);
        edge_modifications_map.update_edge_none(delet.x, delet.y, EdgeType::UNDIRECTED);
    }

    // H = Ne(y) ∩ Ad(x) \ C
    FlatSet H;
    const auto &neighbors_y = get_neighbors(delet.y);
    const auto &adjacent_x = get_adjacent(delet.x);
    std::ranges::set_intersection(neighbors_y, adjacent_x, std::inserter(H, H.begin()));
    // TODO: can be improved by doing deletion with intersection
    for (int z: delet.C) { H.erase(z); }

    // 2. for each h ∈ H:
    //   - orient the (previously undirected) edges between h and y as y → h
    //      [they are all undirected]
    //   - orient any (previously undirected) edges between x and h as x → h
    //      [some are undirected]
    for (int h: H) {
        remove_undirected_edge(delet.y, h);
        add_directed_edge(delet.y, h);
        edge_modifications_map.update_edge_directed(delet.y, h, EdgeType::UNDIRECTED);

        if (has_undirected_edge(delet.x, h)) {
            remove_undirected_edge(delet.x, h);
            add_directed_edge(delet.x, h);
            edge_modifications_map.update_edge_directed(delet.x, h, EdgeType::UNDIRECTED);
        }
    }
    EdgeQueueSet edges_to_check;
    add_adjacent_edges(delet.x, delet.y, edges_to_check);
    for (int h: H) {
        add_adjacent_edges(h, delet.y, edges_to_check);
        add_adjacent_edges(delet.x, h, edges_to_check);
    }
    complete_cpdag(edges_to_check, edge_modifications_map);
    statistics["apply_delete-time"] += measure_time(start_time);
}


/**
 * @brief Check if Meek rule 1 applies:  (z → x - y)  ⇒  (x → y)
 *
 * Condition: Exists z such that
 *  1. z → x
 *  2. z not adjacent to y
 *
 *  @param x The source node
 *  @param y The destination node
 *  @return `true` if Meek rule 1 directs the edge x → y, `false` otherwise.
 */
bool PDAG::is_oriented_by_meek_rule_1(const int x, const int y) const {
    // 1. z → x
    for (int z: parents.at(x)) {
        // 2. z not adjacent to y
        if (!adjacent.at(y).contains(z)) { return true; }
    }
    return false;
}

/**
 * @brief Check if Meek rule 2 applies: (x → z → y) ∧ (x - y)  ⇒  (x → y)
 *
 * Condition: Exists z such that
 *  1. x → z
 *  2. z → y
 *
 *  @param x The source node
 *  @param y The destination node
 *  @return `true` if Meek rule 2 directs the edge x → y, `false` otherwise.
 *
 */
bool PDAG::is_oriented_by_meek_rule_2(const int x, const int y) const {
    // 1. x → z
    for (const int z: children.at(x)) {
        // 2. z → y
        if (children.at(z).contains(y)) { return true; }
    }
    return false;
}

/**
 * @brief Check if Meek rule 3 applies: (x - z → y) ∧ (x - w → y) ∧ (x - y)  ⇒  (x → y)
 *
 * Condition: Exists (z, w) such that
 *  1. z - x and w - x
 *  2. z → y and w → y
 *  3. z ≠ w
 *  4. z, w not adjacent
 *
 *  @param x The source node
 *  @param y The destination node
 *  @return `true` if Meek rule 3 directs the edge x → y, `false` otherwise.
 */
bool PDAG::is_oriented_by_meek_rule_3(int x, int y) const {
    // 1. z - x and w - x
    const auto &neighbors_x = neighbors.at(x);
    const auto &parent_y = parents.at(y);
    FlatSet candidates_z_w;
    // 2. z → y and w → y
    std::ranges::set_intersection(neighbors_x, parent_y,
                                  std::inserter(candidates_z_w, candidates_z_w.begin()));
    for (auto candidate_z: candidates_z_w) {
        for (auto candidate_w: candidates_z_w) {
            // 3. z ≠ w
            if (candidate_z >= candidate_w) { continue; }
            // 4. z, w not adjacent
            auto &adjacent_z = adjacent.at(candidate_z);
            if (adjacent_z.find(candidate_w) == adjacent_z.end()) { return true; }
        }
    }
    return false;
}


/**
 * NOT TESTED; NOT USED
 * @brief Check if Meek rule 4 applies: (w - x - y) ∧ (w → z → y) ∧ (w - y)  ⇒  (x → y)
 *
 * Condition: Exists (z, w) such that
 *  1. w - x and w - y
 *  2. w → z and z → y
 *  3. z, x not adjacent
 */
bool PDAG::is_oriented_by_meek_rule_4(int x, int y) const {
    // 1. w - x and w - y
    const auto &neighbors_x = neighbors.at(x);
    const auto &neighbors_y = neighbors.at(y);
    FlatSet candidates_w;
    std::ranges::set_intersection(neighbors_x, neighbors_y,
                                  std::inserter(candidates_w, candidates_w.begin()));
    for (auto candidate_w: candidates_w) {
        // 2. w → z and z → y
        for (auto candidate_z: children.at(candidate_w)) {
            if (children.at(candidate_z).find(y) != children.at(candidate_z).end()) {
                // 3. z, x not adjacent
                if (adjacent.at(candidate_z).find(x) == adjacent.at(candidate_z).end()) {
                    return true;
                }
            }
        }
    }
    return false;
}

/**
 * @brief Check if (x, y) is part of a v-structure.
 *
 * Conditions: Exists z such that
 *  1. x → y
 *  2. z → y
 *  3. x ≠ z
 *  4. z not adjacent to x.
 *
 * @param x
 * @param y
 * @return `true` if (x, y) is part of a v-structure, `false` otherwise.
 */
bool PDAG::is_part_of_v_structure(int x, int y) const {
    auto &parents_y = parents.at(y);
    //1. x → y
    if (parents_y.find(x) == parents_y.end()) { return false; }
    // 2. z → y
    for (int z: parents_y) {
        // 3. x ≠ z
        if (z == x) { continue; }
        // 4. z not adjacent to x.
        if (adjacent.at(z).find(x) == adjacent.at(z).end()) { return true; }
    }
    return false;
}

/**
 * @brief Complete the PDAG as a CPDAG efficiently.
 *
 * Checks if the edges in `edges_to_check` should be directed or undirected, and update
 * the PDAG accordingly. The modified edges are added to `edge_modifications_map`.
 *
 * Each edge (x, y) in `edges_to_check` is directed if and only if one of the following holds:
 * 1. x → y is part of a v-structure;
 * 2.1 x → y is oriented by Meek rule 1;
 * 2.2 x → y is oriented by Meek rule 2;
 * 2.3. x → y is oriented by Meek rule 3;
 * Otherwise, the edge is undirected.
 *
 * If an edge is updated, the adjacent edges are added to `edges_to_check` to
 * check if they should be updated as well.
 *
 * @param edges_to_check The edges to check
 * @param edge_modifications_map The map of edge modifications
 */
void PDAG::complete_cpdag(EdgeQueueSet &edges_to_check,
                          EdgeModificationsMap &edge_modifications_map) {
    while (!edges_to_check.empty()) {
        Edge edge = edges_to_check.pop();

        if (edge.is_directed()) {
            // Check if the edge is still directed
            int x = edge.get_source();
            int y = edge.get_target();
            if (is_part_of_v_structure(x, y) || is_oriented_by_meek_rule_1(x, y) ||
                is_oriented_by_meek_rule_2(x, y) || is_oriented_by_meek_rule_3(x, y)) {
                // The edge is still directed
                continue;
            }
            // The edge is not directed anymore
            remove_directed_edge(x, y);
            add_undirected_edge(x, y);
            edge_modifications_map.update_edge_undirected(x, y, EdgeType::DIRECTED_TO_Y);
            add_adjacent_edges(x, y, edges_to_check);
        } else {
            // Check if the edge is now directed
            int x = edge.get_x();
            int y = edge.get_y();
            if (is_oriented_by_meek_rule_1(x, y) || is_oriented_by_meek_rule_2(x, y) ||
                is_oriented_by_meek_rule_3(x, y)) {
                // The edge is now directed
            } else if (is_oriented_by_meek_rule_1(y, x) ||
                       is_oriented_by_meek_rule_2(y, x) ||
                       is_oriented_by_meek_rule_3(y, x)) {
                // The edge is now directed
                std::swap(x, y);
            } else {
                // The edge is still undirected
                continue;
            }
            // The edge is now directed
            remove_undirected_edge(x, y);
            add_directed_edge(x, y);
            // the edge was undirected
            edge_modifications_map.update_edge_directed(x, y, EdgeType::UNDIRECTED);
            add_adjacent_edges(x, y, edges_to_check);
        }
    }
    assert(check_is_cpdag());
}

bool PDAG::check_is_cpdag() {
    for (int x: nodes_variables) {
        for (int y: children.at(x)) {
            // check x → y is indeed a directed edge
            if (!(is_part_of_v_structure(x, y) || is_oriented_by_meek_rule_1(x, y) ||
                  is_oriented_by_meek_rule_2(x, y) || is_oriented_by_meek_rule_3(x, y))) {
                return false;
            }
        }
        for (int y: neighbors.at(x)) {
            // check x - y is indeed an undirected edge
            if (is_part_of_v_structure(x, y) || is_oriented_by_meek_rule_1(x, y) ||
                is_oriented_by_meek_rule_2(x, y) || is_oriented_by_meek_rule_3(x, y)) {
                return false;
            }
            if (is_part_of_v_structure(y, x) || is_oriented_by_meek_rule_1(y, x) ||
                is_oriented_by_meek_rule_2(y, x) || is_oriented_by_meek_rule_3(y, x)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Add edges adjacent to (x, y) to the `edgeQueueSet`.
 *
 * @param x The source node
 * @param y The destination node
 * @param edgeQueueSet The set to add the edges to
*/
void PDAG::add_adjacent_edges(int x, int y, EdgeQueueSet &edgeQueueSet) {
    for (int z: children.at(y)) {
        if (x != z) edgeQueueSet.push_directed(y, z);
    }
    for (int z: parents.at(y)) {
        if (x != z) edgeQueueSet.push_directed(z, y);
    }
    for (int z: neighbors.at(y)) {
        if (x != z) edgeQueueSet.push_undirected(y, z);
    }
    for (int z: children.at(x)) {
        if (y != z) edgeQueueSet.push_directed(x, z);
    }
    for (int z: parents.at(x)) {
        if (y != z) edgeQueueSet.push_directed(z, x);
    }
    for (int z: neighbors.at(x)) {
        if (y != z) edgeQueueSet.push_undirected(x, z);
    }
}

PDAG PDAG::get_dag_extension() const {
    PDAG dag_extension = *this;
    PDAG dag_tmp = *this;
    std::set nodes_tmp(nodes_variables.begin(), nodes_variables.end());

    while (!nodes_tmp.empty()) {
        // find a node x that:
        // 1. has no children (children[x] = ∅)
        // 2. For every neighbor y of x, y is adjacent to all the other vertices which are adjacent to x;
        // ∀y ∈ Ne(x) : Ad(x)\{y} ⊆ Ad(y)  i.e. ∀y ∈ Ne(x) : Ad(x) ⊆ Ad(y) ∪ {y}
        int x = -1;

        for (int node: nodes_tmp) {
            if (dag_tmp.get_children(node).empty()) {
                bool is_dag_extension = true;
                for (int neighbor: dag_tmp.get_neighbors(node)) {
                    auto adjacent_neighbor = dag_tmp.get_adjacent(neighbor);
                    adjacent_neighbor.insert(neighbor);
                    if (!is_subset(dag_tmp.get_adjacent(node), adjacent_neighbor)) {
                        is_dag_extension = false;
                        break;
                    }
                }
                if (is_dag_extension) {
                    x = node;
                    break;
                }
            }
        }
        if (x == -1) { throw std::runtime_error("No consistent extension possible"); }
        // Let all the edges which are adjacent to x in dag_tmp be directed toward x in dag_extension
        // node_tmp := node_tmp - x
        // dag_tmp: remove all edges incident to x
        // Have to be very careful with iterators here
        while (!dag_tmp.get_neighbors(x).empty()) {
            int neighbor = *dag_tmp.get_neighbors(x).begin();
            dag_tmp.remove_undirected_edge(neighbor, x);
            dag_extension.remove_undirected_edge(neighbor, x);
            dag_extension.add_directed_edge(neighbor, x);
        }

        while (!dag_tmp.get_parents(x).empty()) {
            int parent = *dag_tmp.get_parents(x).begin();
            dag_tmp.remove_directed_edge(parent, x);
        }
        nodes_tmp.erase(x);
    }
    return dag_extension;
}

std::string PDAG::get_adj_string() const {
    std::string result;
    // first line, each node
    for (int node: nodes_variables) { result += std::to_string(node) + ", "; }
    // remove last ", "
    if (!result.empty()) {
        result.pop_back();
        result.pop_back();
    }
    result += "\n";
    // other line: adjacency matrix (0,1)
    for (int node: nodes_variables) {
        std::string line;
        for (int node2: nodes_variables) {
            if (node == node2) {
                line += "0, ";
            } else if (has_undirected_edge(node, node2) ||
                       has_directed_edge(node, node2)) {
                line += "1, ";
            } else {
                line += "0, ";
            }
        }
        // remove last ", "
        if (!line.empty()) {
            line.pop_back();
            line.pop_back();
        }
        result += line;
        result += "\n";
    }
    return result;
}

PDAG PDAG::from_file(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<int> nodes;

    // Read first line to get nodes
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        int node;
        while (ss >> node) {
            nodes.push_back(node);
            if (ss.peek() == ',') ss.ignore();
        }
    }

    PDAG graph(nodes.size());

    // Read adjacency matrix lines
    int i = 0;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        int j = 0;
        while (std::getline(ss, token, ',')) {
            int value = std::stoi(token);
            if (value == 1 && i != j) { graph.add_directed_edge(nodes[i], nodes[j]); }
            j++;
        }
        i++;
    }
    file.close();
    return graph;
}

/**
 * @brief Compute the SHD score between two PDAGs.
 *
 * The SHD score is the number of edges that are different between the two PDAGs.
 * If `allow_directed_in_other` is `true`, the score is reduced by 1 for each undirected edge
 * in the first PDAG that is directed in the second PDAG.
 *
 * @param other The other PDAG
 * @param allow_directed_in_other If directed edges in the second PDAG should be allowed
 * @return The SHD score
 */
int PDAG::shd(const PDAG &other, bool allow_directed_in_other) const {
    int shd = 0;
    for (int node: nodes_variables) {
        for (int node2: nodes_variables) {
            if (node >= node2) { continue; }
            if ((has_directed_edge(node, node2) &&
                 !other.has_directed_edge(node, node2)) ||
                (has_directed_edge(node2, node) &&
                 !other.has_directed_edge(node2, node)) ||
                (has_undirected_edge(node, node2) &&
                 !other.has_undirected_edge(node, node2)) ||
                (adjacent.at(node).find(node2) == adjacent.at(node).end() &&
                 other.adjacent.at(node).find(node2) != other.adjacent.at(node).end())) {
                shd++;
            }
            if (allow_directed_in_other && has_undirected_edge(node, node2) &&
                (other.has_directed_edge(node, node2) ||
                 other.has_directed_edge(node2, node))) {
                shd--;
            }
        }
    }
    return shd;
}

std::ostream &operator<<(std::ostream &os, const PDAG &obj) {
    os << "PDAG: undirected edges = {";
    for (auto [fst, snd]: obj.get_undirected_edges()) {
        os << "(" << fst << "-" << snd << "), ";
    }
    os << "}, directed edges = {";
    for (auto [fst, snd]: obj.get_directed_edges()) {
        os << "(" << fst << "→" << snd << "), ";
    }
    os << "}";
    return os;
}

bool PDAG::operator==(const PDAG &other) const {
    if (number_of_directed_edges != other.number_of_directed_edges ||
        number_of_undirected_edges != other.number_of_undirected_edges) {
        return false;
    }
    if (nodes_variables.size() != other.nodes_variables.size()) { return false; }
    for (int node: nodes_variables) {
        if (children.at(node) != other.children.at(node) ||
            parents.at(node) != other.parents.at(node) ||
            neighbors.at(node) != other.neighbors.at(node)) {
            return false;
        }
    }
    return true;
}
