//
// Created by Achille Nazaret on 11/3/23.
//
#include "PDAG.h"
#include "set_ops.h"
#include <fstream>
#include <sstream>

using namespace std::chrono;

PDAG::PDAG(int num_variables, int num_interventions)
    : num_variables(num_variables), num_interventions(num_interventions),
      _block_semi_directed_path_queue(num_variables) {
    for (int i = 0; i < num_variables; i++) {
        nodes_variables.push_back(i);
        nodes_all.push_back(i);
        children.emplace_back();
        parents.emplace_back();
        neighbors.emplace_back();
        adjacent.emplace_back();
        adjacent_reachable.emplace_back();
        node_version.push_back(0);
    }
    for (int i = 0; i < num_interventions; i++) {
        nodes_interventions.push_back(i + num_variables);
        nodes_all.push_back(i + num_variables);
        children.emplace_back();
        parents.emplace_back();
        neighbors.emplace_back();
        adjacent.emplace_back();
        adjacent_reachable.emplace_back();
        node_version.push_back(0);
    }

    _block_semi_directed_path_visited.resize(num_variables);
    _block_semi_directed_path_blocked.resize(num_variables);
}

int PDAG::get_number_of_edges() const { return number_of_directed_edges + number_of_undirected_edges; }

int PDAG::get_node_version(int node) const { return node_version.at(node); }

const std::vector<int> &PDAG::get_nodes_variables() const { return nodes_variables; }

std::vector<std::pair<int, int>> PDAG::get_directed_edges() const {
    std::vector<std::pair<int, int>> edges;
    for (const auto &node: nodes_variables) {
        for (int child: children.at(node)) { edges.emplace_back(node, child); }
    }
    return edges;
}

std::vector<std::pair<int, int>> PDAG::get_directed_edges_interventions() const {
    std::vector<std::pair<int, int>> edges;
    for (const auto &node: nodes_interventions) {
        for (int child: children.at(node)) { edges.emplace_back(node, child); }
    }
    return edges;
}

std::vector<std::pair<int, int>> PDAG::get_undirected_edges() const {
    std::vector<std::pair<int, int>> edges;
    // interventions cannot have undirected edge
    for (const auto &node: nodes_variables) {
        for (int neighbor: neighbors.at(node)) {
            if (node < neighbor) { edges.emplace_back(node, neighbor); }
        }
    }
    return edges;
}

const FlatSet &PDAG::get_parents(int node) const { return parents[node]; }

const FlatSet &PDAG::get_children(int node) const { return children[node]; }

const FlatSet &PDAG::get_neighbors(int node) const { return neighbors[node]; }

const FlatSet &PDAG::get_adjacent(int node) const { return adjacent[node]; }

const FlatSet &PDAG::get_adjacent_reachable(int node) const { return adjacent_reachable.at(node); }

FlatSet PDAG::get_neighbors_adjacent(int node_y, int node_x) const {
    FlatSet result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);
    std::set_intersection(neighbors_set.begin(), neighbors_set.end(), adjacent_set.begin(),
                          adjacent_set.end(), std::inserter(result, result.begin()));
    return result;
}

FlatSet PDAG::get_neighbors_not_adjacent(int node_y, int node_x) const {
    FlatSet result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);
    std::set_difference(neighbors_set.begin(), neighbors_set.end(), adjacent_set.begin(), adjacent_set.end(),
                        std::inserter(result, result.begin()));
    return result;
}

bool PDAG::is_clique(const FlatSet &nodes_subset) const {
    for (int node1: nodes_subset) {
        const FlatSet &adjacent_set = get_adjacent(node1);
        for (int node2: nodes_subset) {
            if (node1 != node2 && adjacent_set.find(node2) == adjacent_set.end()) { return false; }
        }
    }
    return true;
}

bool PDAG::has_directed_edge(int x, int y) const {
    auto &children_set = children.at(x);
    return children_set.find(y) != children_set.end();
}

bool PDAG::has_undirected_edge(int x, int y) const {
    auto &neighbors_set = neighbors.at(x);
    return neighbors_set.find(y) != neighbors_set.end();
}

void PDAG::remove_directed_edge(int x, int y) {
    children[x].erase(y);
    parents[y].erase(x);
    adjacent[x].erase(y);
    adjacent[y].erase(x);
    adjacent_reachable[x].erase(y);
    number_of_directed_edges--;
    node_version[x]++;
    node_version[y]++;
    graph_version++;
    modifications_history.emplace_back(PDAGModification::REMOVE_DIRECTED_EDGE, x, y);
}

void PDAG::remove_undirected_edge(int x, int y) {
    neighbors[x].erase(y);
    neighbors[y].erase(x);
    adjacent[x].erase(y);
    adjacent[y].erase(x);
    adjacent_reachable[x].erase(y);
    adjacent_reachable[y].erase(x);
    number_of_undirected_edges--;
    node_version[x]++;
    node_version[y]++;
    graph_version++;
    modifications_history.emplace_back(PDAGModification::REMOVE_UNDIRECTED_EDGE, x, y);
}

void PDAG::add_directed_edge(int x, int y) {
    children[x].insert(y);
    parents[y].insert(x);
    adjacent[x].insert(y);
    adjacent[y].insert(x);
    adjacent_reachable[x].insert(y);
    number_of_directed_edges++;
    node_version[x]++;
    node_version[y]++;
    graph_version++;
    modifications_history.emplace_back(PDAGModification::ADD_DIRECTED_EDGE, x, y);
}

void PDAG::add_undirected_edge(int x, int y) {
    neighbors[x].insert(y);
    neighbors[y].insert(x);
    adjacent[x].insert(y);
    adjacent[y].insert(x);
    adjacent_reachable[x].insert(y);
    adjacent_reachable[y].insert(x);
    number_of_undirected_edges++;
    node_version[x]++;
    node_version[y]++;
    graph_version++;
    modifications_history.emplace_back(PDAGModification::ADD_UNDIRECTED_EDGE, x, y);
}

/**
 * Check if blocked_nodes block all semi-directed paths from src to dst.
 *
 * @param src The source node
 * @param dst The destination node
 * @param blocked_nodes The set of nodes that are blocked
 * @return `true` if blocked_nodes block all semi-directed paths from src to dst, `false` otherwise.
 */
bool PDAG::block_semi_directed_paths(int src, int dst, const FlatSet &blocked_nodes,
                                     bool ignore_direct_edge) {
    statistics["call- block_semi_directed_paths"] += 1;
    if (src == dst) { return false; }
    // BFS search from y to x, using adjacent_reachable edges, avoiding blocked nodes
    auto &visited = _block_semi_directed_path_visited;
    std::fill(visited.begin(), visited.end(), 0);
    auto &blocked = _block_semi_directed_path_blocked;
    std::fill(blocked.begin(), blocked.end(), 0);

    for (int n: blocked_nodes) { blocked[n] = 1; }

    visited[src] = 1;

    auto &queue = _block_semi_directed_path_queue;
    queue.clear();
    queue.push_back(src);

    while (!queue.empty()) {
        int node = queue.pop_front();
        auto &reachable = get_adjacent_reachable(node);

        for (int n: reachable) {
            if (visited[n] || blocked[n] || (node == src && n == dst && ignore_direct_edge)) { continue; }
            if (n == dst) { return false; }
            queue.push_back(n);
            visited[n] = 1;
        }
    }
    return true;
}


/**
 * Verify if the insert is valid.
 *
 * Insert(x, y, T) is valid if and only if: [in approximate order of complexity]
 *  1. x and y are not adjacent
 *  2. T is a subset of Ne(y) \ Ad(x)
 *  3. The score has not changed, i.e.  [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y) has not changed
 *  4. [Ne(y) ∩ Ad(x)] ∪ T is a clique (even if 3. hold, the clique might not be valid)
 *  5. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x
 *
 * @param insert
 * @return ??
 */
bool PDAG::is_insert_valid(const Insert &insert, bool reverse) {
    statistics["call- is_insert_valid"] += 1;
    int x = insert.x;
    int y = insert.y;
    auto &T = insert.T;

    // add version thingy

    auto &adjacent_x = adjacent.at(x);
    if (!reverse) {
        // 1. x and y are not adjacent
        if (adjacent_x.find(y) != adjacent_x.end()) {
            statistics["is_insert_valid: false 1a"] += 1;
            return false;
        }
    } else {
        // 1. x ← y
        auto &parents_x = parents.at(x);
        if (parents_x.find(y) == parents_x.end()) {
            statistics["is_insert_valid: false 1b"] += 1;
            return false;
        }
    }

    // 2. T ⊆ Ne(y) \ Ad(x)
    // <=> T ⊆ Ne(y) and T does not intersect Ad(x)
    auto &neighbors_y = neighbors.at(y);
    if (!is_subset(T, neighbors_y) || have_overlap(T, adjacent_x)) {
        statistics["is_insert_valid: false 2"] += 1;
        return false;
    }

    // 3. [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y) has not changed
    // <=> insert.effective_parents == [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    FlatSet ne_y_ad_x_T;
    // intersect effective_parents with Ad(x)
    std::set_intersection(neighbors_y.begin(), neighbors_y.end(), adjacent_x.begin(), adjacent_x.end(),
                          std::inserter(ne_y_ad_x_T, ne_y_ad_x_T.begin()));
    // todo: can easily be improved by doing the insertion manually with intersection
    ne_y_ad_x_T.insert(T.begin(), T.end());
    if (!equal_union(insert.effective_parents, ne_y_ad_x_T, parents.at(y))) {
        statistics["is_insert_valid: false 3"] += 1;
        return false;
    }

    // 4. [Ne(y) ∩ Ad(x)] ∪ T is a clique
    if (!is_clique(ne_y_ad_x_T)) {
        statistics["is_insert_valid: false 4"] += 1;
        return false;
    }

    // 5. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x
    bool ignore_direct_edge = reverse;
    auto start_time = high_resolution_clock::now();
    if (!block_semi_directed_paths(y, x, ne_y_ad_x_T, ignore_direct_edge)) {
        statistics["is_insert_valid: false 5"] += 1;
        statistics["time- block_semi_directed_paths: false"] +=
                duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();
        return false;
    }
    statistics["time- block_semi_directed_paths: true"] +=
            duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();

    return true;
}

bool PDAG::is_reverse_valid(const Reverse &reverse) {
    // TODO: wrong, need to check if the score to x has changed
    // is Pa(x) unchanged
    int x = reverse.insert.x;
    if (get_parents(x) != reverse.parents_x) { return false; }
    return is_insert_valid(reverse.insert, true);
}

bool PDAG::is_delete_valid(const Delete &delet) const {
    // 1. x and y are neighbors or x is a parent of y [aka y is adjacent_reachable from x]
    int x = delet.x;
    int y = delet.y;
    if (adjacent_reachable.at(x).find(y) == adjacent_reachable.at(x).end()) { return false; }

    // 2. O is a subset of Ne(y) ∩ Ad(x)
    // <=> O ⊆ Ne(y) and O ⊆ Ad(x)
    auto &neighbors_y = neighbors.at(y);
    auto &adjacent_x = adjacent.at(x);
    if (!is_subset(delet.O, neighbors_y) || !is_subset(delet.O, adjacent_x)) { return false; }

    // 3. O and Pa(y) are unchanged (technically only need O ∪ Pa(y) to be
    // unchanged)
    //    // <=> delet.effective_parents == O ∪ Pa(y) ∪ {x}
    if (!equal_union_with_singleton(delet.effective_parents, delet.O, parents.at(y), x)) { return false; }

    // 4. O is a clique
    if (!is_clique(delet.O)) { return false; }

    return true;
}

/**
 * Apply the insert to the PDAG, and adapt the graph to remain a CPDAG.
 *
 * Insert(x, y, T) does:
 *  1. insert the directed edge x → y
 *  2. for each t ∈ T: orient the (previously undirected) edge between t and y as t → y
 *
 * @param insert
 */
void PDAG::apply_insert(const Insert &insert, EdgeModificationsMap &edge_modifications_map) {
    int x = insert.x;
    int y = insert.y;
    auto &T = insert.T;


    // Start: The graph is a CPDAG
    EdgeQueueSet edges_to_check;

    // Case 1: T is not empty;
    //  then x → y and each t → y will form a v-structure [remember t ∉ Ad(x)];
    //  so the edges are directed (and we do not need to check them later).
    if (!T.empty()) {
        // 1. insert the directed edge x → y
        add_directed_edge(x, y);
        edge_modifications_map.update_edge_directed(x, y, EdgeType::NONE);

        // 2. for each t ∈ T: orient the (previously undirected) edge between t and y as t → y
        for (int t: T) {
            remove_undirected_edge(t, y);
            add_directed_edge(t, y);
            edge_modifications_map.update_edge_directed(t, y, EdgeType::UNDIRECTED);
        }

        // only now:
        add_edges_around(x, y, edges_to_check, true);
        for (int t: T) { add_edges_around(t, y, edges_to_check, true); }
    } else {
        // Case 2: T is empty;
        //  then x → y is part of a v-structure if and only if Pa(y) \ Ad(x) ≠ ∅.
        //  <==> ¬(Pa(y) ⊆ Ad(x))
        if (!is_subset(parents.at(y), adjacent.at(x))) {
            // actually we want to know the content of Pa(y) \ Ad(x)
            // 1. insert the directed edge x → y
            add_directed_edge(x, y);
            edges_to_check.push_directed(x, y);
            add_edges_around(x, y, edges_to_check, true);
            edge_modifications_map.update_edge_directed(x, y, EdgeType::NONE);
        } else {
            // 2. insert the undirected edge x - y
            add_undirected_edge(x, y);
            edges_to_check.push_undirected(x, y);
            add_edges_around(x, y, edges_to_check, false);
            edge_modifications_map.update_edge_undirected(x, y, EdgeType::NONE);
        }
    }

    // now need to decide which edges to check,

    // i have added edge (x,y)
    // i need to check the adjacent edges that might be affected by this update
    // easy answer: all edges adjacent to x or y
    // but i can do better
    // If x → y, then don't need to update the edges going to x

    maintain_cpdag(edges_to_check, edge_modifications_map);

    // End: The graph is a CPDAG
    // print changed edges
    //    std::cout << "changed_edges.size() = " << changed_edges.size() << ": ";
    //    for (auto edge: changed_edges) {
    //        std::cout << edge.x << (edge.type == EdgeType::DIRECTED ? " → " : " - ") << edge.y << ", ";
    //    }
    //    std::cout << std::endl;
}

void PDAG::apply_reverse(const Reverse &reverse, EdgeModificationsMap &edge_modifications_map) {
    // Start: The graph is a CPDAG

    // todo: if we change the tracking of the edges, we'll have to add y → x here
    // 1. remove the directed edge y → x
    int x = reverse.insert.x;
    int y = reverse.insert.y;
    remove_directed_edge(y, x);
    edge_modifications_map.update_edge_none(x, y, EdgeType::DIRECTED_TO_X);
    apply_insert(reverse.insert, edge_modifications_map);
}

void PDAG::apply_delete(const Delete &aDelete, EdgeModificationsMap &edge_modifications_map) {
    // Start: The graph is a CPDAG
    if (has_directed_edge(aDelete.x, aDelete.y)) {
        // 1. remove the directed edge x → y
        remove_directed_edge(aDelete.x, aDelete.y);
        edge_modifications_map.update_edge_none(aDelete.x, aDelete.y, EdgeType::DIRECTED_TO_Y);
    } else {
        // 1. remove the undirected edge x - y
        remove_undirected_edge(aDelete.x, aDelete.y);
        edge_modifications_map.update_edge_none(aDelete.x, aDelete.y, EdgeType::UNDIRECTED);
    }

    // H = Ne(y) ∩ Ad(x) \ O
    FlatSet H;
    const auto &neighbors_y = get_neighbors(aDelete.y);
    const auto &adjacent_x = get_adjacent(aDelete.x);
    std::set_intersection(neighbors_y.begin(), neighbors_y.end(), adjacent_x.begin(), adjacent_x.end(),
                          std::inserter(H, H.begin()));
    // TODO: can easily be improved by doing deletion with intersection
    for (int z: aDelete.O) { H.erase(z); }


    // 2. for each h ∈ H:
    //   - orient the (previously undirected) edge between h and y as y → h [they are all undirected]
    //   - orient any (previously undirected) edge between x and h as x → h [some might be undirected]
    for (int h: H) {
        remove_undirected_edge(aDelete.y, h);
        add_directed_edge(aDelete.y, h);
        edge_modifications_map.update_edge_directed(aDelete.y, h, EdgeType::UNDIRECTED);

        if (has_undirected_edge(aDelete.x, h)) {
            remove_undirected_edge(aDelete.x, h);
            add_directed_edge(aDelete.x, h);
            edge_modifications_map.update_edge_directed(aDelete.x, h, EdgeType::UNDIRECTED);
        }
    }
    EdgeQueueSet edges_to_check;

    add_edges_around(aDelete.x, aDelete.y, edges_to_check, false);
    for (int h: H) {
        add_edges_around(h, aDelete.y, edges_to_check, true);
        add_edges_around(aDelete.x, h, edges_to_check, true);// shouldnt make a difference
    }
    maintain_cpdag(edges_to_check, edge_modifications_map);
}


void PDAG::add_edges_around(int x, int y, EdgeQueueSet &edgeQueueSet, bool is_directed) {
    for (int z: children.at(y)) {
        if (x != z) edgeQueueSet.push_directed(y, z);
    }
    for (int z: parents.at(y)) {
        if (x != z) edgeQueueSet.push_directed(z, y);
    }
    for (int z: children.at(x)) {
        if (y != z) edgeQueueSet.push_directed(x, z);
    }
    // if is_directed, we do not need to add any (z, x) edges for z in parents.at(x)
    if (!is_directed) {
        for (int z: parents.at(x)) {
            if (y != z) edgeQueueSet.push_directed(z, x);
        }
    }
    for (int z: neighbors.at(x)) {
        if (y != z) edgeQueueSet.push_undirected(x, z);
    }
    for (int z: neighbors.at(y)) {
        if (x != z) edgeQueueSet.push_undirected(y, z);
    }
}

/**
 * Meek rule 1:  (z → x - y)  ⇒  (x → y)
 *
 * Assume that we have (x - y)
 * Condition: Exists z such that
 *  1. z → x
 *  2. z not adjacent to y
 */
bool PDAG::is_oriented_by_meek_rule_1(int x, int y) const {
    // 1. z → x
    for (int z: parents.at(x)) {
        // 2. z not adjacent to y
        if (adjacent.at(z).find(y) == adjacent.at(z).end()) { return true; }
    }
    return false;
}

/**
 * Meek rule 2: (x → z → y) ∧ (x - y)  ⇒  (x → y)
 *
 * Assume that we have (x - y)
 * Condition: Exists z such that
 *  1. x → z
 *  2. z → y
 */
bool PDAG::is_oriented_by_meek_rule_2(int x, int y) const {
    // 1. x → z
    for (int z: children.at(x)) {
        // 2. z → y
        if (children.at(z).find(y) != children.at(z).end()) { return true; }
    }
    return false;
}

/**
 * Meek rule 3: (x - z → y) ∧ (x - w → y) ∧ (x - y)  ⇒  (x → y)
 *
 * Assume that we have (x - y)
 * Condition: Exists (z, w) such that
 *  1. z - x and w - x
 *  2. z → y and w → y
 *  3. z ≠ w
 *  4. z, w not adjacent
 */
bool PDAG::is_oriented_by_meek_rule_3(int x, int y) const {
    // 1. z - x and w - x
    const auto &neighbors_x = neighbors.at(x);
    const auto &parent_y = parents.at(y);
    FlatSet candidates_z_w;
    // 2. z → y and w → y
    std::set_intersection(neighbors_x.begin(), neighbors_x.end(), parent_y.begin(), parent_y.end(),
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
 * NOT TESTED
 * Meek rule 4: (w - x - y) ∧ (w → z → y) ∧ (w - y)  ⇒  (x → y)
 *
 * Assume that we have (x - y)
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
    std::set_intersection(neighbors_x.begin(), neighbors_x.end(), neighbors_y.begin(), neighbors_y.end(),
                          std::inserter(candidates_w, candidates_w.begin()));
    for (auto candidate_w: candidates_w) {
        // 2. w → z and z → y
        for (auto candidate_z: children.at(candidate_w)) {
            if (children.at(candidate_z).find(y) != children.at(candidate_z).end()) {
                // 3. z, x not adjacent
                if (adjacent.at(candidate_z).find(x) == adjacent.at(candidate_z).end()) { return true; }
            }
        }
    }
    return false;
}

/**
 * Check if (x, y) is part of a v-structure.
 *
 * Condition: Exists z such that:
 *  1. x → y
 *  2. z → y
 *  3. x ≠ z
 *  4. z not adjacent to x.
 *
 * @param x
 * @param y
 * @return
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
 * Check the types of a list of edges.
 *
 * Each edge (x, y) in edges_to_check is directed if and only if one of the following holds:
 * 1. x → y is part of a v-structure;
 * 2.1 x → y is oriented by Meek rule 1;
 * 2.2 x → y is oriented by Meek rule 2;
 * 2.3. x → y is oriented by Meek rule 3;
 *
 * Otherwise, the edge is undirected.
 *
 * @param insert
 */
void PDAG::maintain_cpdag(EdgeQueueSet &edges_to_check, EdgeModificationsMap &edge_modifications_map) {
    // need to keep track of the modifications made to the graph

    while (!edges_to_check.empty()) {
        Edge edge = edges_to_check.pop();
        int x = edge.x;
        int y = edge.y;
        bool new_is_directed;

        if (edge.type == EdgeType::DIRECTED_TO_Y) {
            // Check if the edge is still directed
            if (is_part_of_v_structure(x, y) || is_oriented_by_meek_rule_1(x, y) ||
                is_oriented_by_meek_rule_2(x, y) || is_oriented_by_meek_rule_3(x, y) ||
                node_is_intervention(x)) {
                // The edge is still directed
                continue;
            }
            // The edge is not directed anymore
            remove_directed_edge(x, y);
            add_undirected_edge(x, y);
            edge_modifications_map.update_edge_undirected(x, y, edge.type);
            new_is_directed = false;
        } else {
            assert(edge.type != EdgeType::DIRECTED_TO_X && edge.type != EdgeType::NONE);
            // Check if the edge is now directed
            if (is_oriented_by_meek_rule_1(x, y) || is_oriented_by_meek_rule_2(x, y) ||
                is_oriented_by_meek_rule_3(x, y) || node_is_intervention(x)) {

                // The edge is now directed
            } else if (is_oriented_by_meek_rule_1(y, x) || is_oriented_by_meek_rule_2(y, x) ||
                       is_oriented_by_meek_rule_3(y, x) || node_is_intervention(y)) {
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
            assert(edge.type == EdgeType::UNDIRECTED);
            edge_modifications_map.update_edge_directed(x, y, edge.type);
            new_is_directed = true;
        }
        // If we reach this point, the edge has been modified
        // We need to check the adjacent edges that might be affected by this update
        add_edges_around(x, y, edges_to_check, new_is_directed);
    }
}

PDAG PDAG::get_dag_extension() const {
    PDAG dag_extension = *this;
    PDAG dag_tmp = *this;
    // todo: check only using nodes_variables is ok
    std::set<int> nodes_tmp(nodes_variables.begin(), nodes_variables.end());

    while (nodes_tmp.size() > 0) {
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
        if (x == -1) {
            // raise error, no consistent extension possible
            std::cout << "no consistent extension possible" << std::endl;
            // todo: raise error
            break;
        }
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
    std::string result = "";
    // first line, each node
    for (int node: nodes_variables) { result += std::to_string(node) + ", "; }
    // remove last ", "
    if (!result.empty()) {
        result.pop_back();
        result.pop_back();
    }
    result += "\n";
    std::string line;
    // other line: adjacency matrix (0,1)
    for (int node: nodes_variables) {
        line = "";
        for (int node2: nodes_variables) {
            if (node == node2) {
                line += "0, ";
            } else if (has_undirected_edge(node, node2) || has_directed_edge(node, node2)) {
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

    PDAG graph(nodes.size(), 0);

    // Read adjacency matrix lines
    int i = 0;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        int j = 0;
        while (std::getline(ss, token, ',')) {
            int value = std::stoi(token);
            if (value == 1 && i != j) {
                graph.add_directed_edge(nodes[i], nodes[j]);
                std::cout << "add_directed_edge(" << nodes[i] << ", " << nodes[j] << ")" << std::endl;
            }
            j++;
        }
        i++;
        std::cout << "i = " << i << std::endl;
    }
    file.close();
    return graph;
}

int PDAG::shd(const PDAG &other, bool allow_directed_in_other) const {
    int shd = 0;
    for (int node: nodes_variables) {
        for (int node2: nodes_variables) {
            if (node >= node2) { continue; }
            if ((has_directed_edge(node, node2) && !other.has_directed_edge(node, node2)) ||
                (has_directed_edge(node2, node) && !other.has_directed_edge(node2, node)) ||
                (has_undirected_edge(node, node2) && !other.has_undirected_edge(node, node2)) ||
                (adjacent.at(node).find(node2) == adjacent.at(node).end() &&
                 other.adjacent.at(node).find(node2) != other.adjacent.at(node).end())) {
                shd++;
            }
            if (allow_directed_in_other && has_undirected_edge(node, node2) &&
                (other.has_directed_edge(node, node2) || other.has_directed_edge(node2, node))) {
                shd--;
            }
        }
    }
    return shd;
}

std::ostream &operator<<(std::ostream &os, const PDAG &obj) {
    os << "PDAG: interventions edges = {";
    for (auto edge: obj.get_directed_edges_interventions()) {
        os << "(I" << edge.first - obj.num_variables << "-->" << edge.second << "), ";
    }
    os << "}, directed edges = {";
    for (auto edge: obj.get_directed_edges()) { os << "(" << edge.first << "-->" << edge.second << "), "; }
    os << "}, undirected edges = {";
    for (auto edge: obj.get_undirected_edges()) { os << "(" << edge.first << "---" << edge.second << "), "; }
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
        if (children.at(node) != other.children.at(node) || parents.at(node) != other.parents.at(node) ||
            neighbors.at(node) != other.neighbors.at(node)) {
            return false;
        }
    }
    return true;
}