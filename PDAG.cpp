//
// Created by Achille Nazaret on 11/3/23.
//
#include "PDAG.h"
#include "set_ops.h"

PDAG::PDAG(int num_nodes) : num_nodes(num_nodes) {
    for (int i = 0; i < num_nodes; i++) {
        nodes.push_back(i);
        children.emplace_back();
        parents.emplace_back();
        neighbors.emplace_back();
        adjacent.emplace_back();
        adjacent_reachable.emplace_back();
        node_version.push_back(0);
    }
}

int PDAG::get_number_of_edges() const { return number_of_directed_edges + number_of_undirected_edges; }

int PDAG::get_node_version(int node) const { return node_version.at(node); }

const std::vector<int> &PDAG::get_nodes() const { return nodes; }

// TODO: Find where this is used and potentially change it to std::vector<std::pair<int, int>>.
std::set<std::pair<int, int>> PDAG::get_directed_edges() const {
    std::set<std::pair<int, int>> edges;
    for (const auto &node: nodes) {
        for (int child: children.at(node)) { edges.insert({node, child}); }
    }
    return edges;
}

// TODO: Find where this is used and potentially change it to std::vector<std::pair<int, int>>.
std::set<std::pair<int, int>> PDAG::get_undirected_edges() const {
    std::set<std::pair<int, int>> edges;
    for (const auto &node: nodes) {
        for (int neighbor: neighbors.at(node)) {
            if (node < neighbor) { edges.insert({node, neighbor}); }
        }
    }
    return edges;
}

const std::set<int> &PDAG::get_parents(int node) const { return parents[node]; }

const std::set<int> &PDAG::get_children(int node) const { return children[node]; }

const std::set<int> &PDAG::get_neighbors(int node) const { return neighbors[node]; }

const std::set<int> &PDAG::get_adjacent(int node) const { return adjacent[node]; }

const std::set<int> PDAG::get_adjacent_reachable(int node) const { return adjacent_reachable.at(node); }

std::set<int> PDAG::get_neighbors_adjacent(int node_y, int node_x) const {
    std::set<int> result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);
    std::set_intersection(neighbors_set.begin(), neighbors_set.end(), adjacent_set.begin(),
                          adjacent_set.end(), std::inserter(result, result.begin()));
    return result;
}

std::set<int> PDAG::get_neighbors_not_adjacent(int node_y, int node_x) const {
    std::set<int> result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);
    std::set_difference(neighbors_set.begin(), neighbors_set.end(), adjacent_set.begin(), adjacent_set.end(),
                        std::inserter(result, result.begin()));
    return result;
}

bool PDAG::is_clique(const std::set<int> &nodes_subset) const {
    for (int node1: nodes_subset) {
        const std::set<int> &adjacent_set = get_adjacent(node1);
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
bool PDAG::block_semi_directed_paths(int src, int dst, const std::set<int> &blocked_nodes) const {
    if (src == dst) { return false; }
    // BFS search from y to x, using adjacent_reachable edges, avoiding blocked nodes
    std::set<int> visited;
    visited.insert(src);

    std::queue<int> queue;
    queue.push(src);

    while (!queue.empty()) {
        int node = queue.front();
        queue.pop();
        // it would be possible to not do a copy and not do find later (probably a minor optimization)
        auto reachable = get_adjacent_reachable(node);

        if (reachable.find(dst) != reachable.end()) { return false; }
        // remove blocked nodes and already visited nodes
        for (auto it = reachable.begin(); it != reachable.end();) {
            if (blocked_nodes.find(*it) != blocked_nodes.end() || visited.find(*it) != visited.end()) {
                it = reachable.erase(it);
            } else {
                ++it;
            }
        }
        for (int n: reachable) {
            queue.push(n);
            visited.insert(n);
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
bool PDAG::is_insert_valid(const Insert &insert) const {
    int x = insert.x;
    int y = insert.y;
    auto &T = insert.T;

    // add version thingy

    // 1. x and y are not adjacent
    auto &adjacent_x = adjacent.at(x);
    if (adjacent_x.find(y) != adjacent_x.end()) { return false; }

    // 2. T ⊆ Ne(y) \ Ad(x)
    // <=> T ⊆ Ne(y) and T does not intersect Ad(x)
    auto &neighbors_y = neighbors.at(y);
    if (!is_subset(T, neighbors_y) || have_overlap(T, adjacent_x)) { return false; }

    // 3. [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y) has not changed
    // <=> insert.effective_parents == [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    auto ne_y_ad_x_T = neighbors_y;
    // intersect effective_parents with Ad(x)
    intersect_in_place(ne_y_ad_x_T, adjacent_x);
    ne_y_ad_x_T.insert(T.begin(), T.end());
    if (!equal_union(insert.effective_parents, ne_y_ad_x_T, parents.at(y))) { return false; }

    // 4. [Ne(y) ∩ Ad(x)] ∪ T is a clique
    if (!is_clique(ne_y_ad_x_T)) { return false; }

    // 5. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x
    if (!block_semi_directed_paths(y, x, ne_y_ad_x_T)) { return false; }

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
void PDAG::apply_insert(const Insert &insert) {
    int x = insert.x;
    int y = insert.y;
    auto &T = insert.T;


    // Start: The graph is a CPDAG
    EdgeQueueSet edges_to_check;
    std::set<Edge> changed_edges;

    // Case 1: T is not empty;
    //  then x → y and each t → y will form a v-structure [remember t ∉ Ad(x)];
    //  so the edges are directed (and we do not need to check them later).
    if (!T.empty()) {
        // 1. insert the directed edge x → y
        add_directed_edge(x, y);
        add_edges_around(x, y, edges_to_check, true);
        changed_edges.insert({x, y, EdgeType::DIRECTED});

        // 2. for each t ∈ T: orient the (previously undirected) edge between t and y as t → y
        for (int t: T) {
            remove_undirected_edge(t, y);
            add_directed_edge(t, y);
            add_edges_around(t, y, edges_to_check, true);
            changed_edges.insert({t, y, EdgeType::DIRECTED});
        }
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
            changed_edges.insert({x, y, EdgeType::DIRECTED});
        } else {
            // 2. insert the undirected edge x - y
            add_undirected_edge(x, y);
            edges_to_check.push_undirected(x, y);
            add_edges_around(x, y, edges_to_check, false);
            changed_edges.insert({x, y, EdgeType::UNDIRECTED});
        }
    }

    // now need to decide which edges to check,

    // Transport to new contexts
    // Transport to new perturbations
    // but trained or not trained?
    // cluster genes together
    // where does the large width regime comes into play?
    // what is the difference between the two regimes?
    // what is the structure of the neural network? how do you get a causal graph Z


    // i have added edge (x,y)
    // i need to check the adjacent edges that might be affected by this update
    // easy answer: all edges adjacent to x or y
    // but i can do better
    // If x → y, then don't need to update the edges going to x

    maintain_cpdag(edges_to_check, changed_edges);

    // End: The graph is a CPDAG
    // print changed edges
    std::cout << "changed_edges.size() = " << changed_edges.size() << std::endl;
    for (auto edge: changed_edges) {
        std::cout << edge.x << (edge.type == EdgeType::DIRECTED ? " → " : " - ") << edge.y << ", ";
    }
    std::cout << std::endl;
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
    auto candidates_z_w = neighbors.at(x);
    // 2. z → y and w → y
    intersect_in_place(candidates_z_w, parents.at(y));
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
void PDAG::maintain_cpdag(EdgeQueueSet &edges_to_check, std::set<Edge> &changed_edges) {
    // need to keep track of the modifications made to the graph

    while (!edges_to_check.empty()) {
        Edge edge = edges_to_check.pop();
        int x = edge.x;
        int y = edge.y;

        if (edge.type == EdgeType::DIRECTED) {
            // Check if the edge is still directed
            if (is_part_of_v_structure(x, y) || is_oriented_by_meek_rule_1(x, y) ||
                is_oriented_by_meek_rule_2(x, y) || is_oriented_by_meek_rule_3(x, y)) {
                // The edge is still directed
                continue;
            }
            // The edge is not directed anymore
            remove_directed_edge(x, y);
            add_undirected_edge(x, y);
        } else {
            // Check if the edge is now directed
            if (is_oriented_by_meek_rule_1(x, y) || is_oriented_by_meek_rule_2(x, y) ||
                is_oriented_by_meek_rule_3(x, y)) {
                // The edge is now directed
            } else if (is_oriented_by_meek_rule_1(y, x) || is_oriented_by_meek_rule_2(y, x) ||
                       is_oriented_by_meek_rule_3(y, x)) {
                // The edge is now directed
                std::swap(x, y);
            } else {
                // The edge is still undirected
                continue;
            }
            remove_undirected_edge(x, y);
            add_directed_edge(x, y);
        }

        // remove the edge with its old type if it was in the list of changed edges
        changed_edges.erase(
                edge);// TODO: if the edge is undirected, we need to check if (v,u) is in changed_edges
        // Add the edge to the list of changed edges, with its new type
        EdgeType new_type = edge.type == EdgeType::UNDIRECTED ? EdgeType::DIRECTED : EdgeType::UNDIRECTED;
        changed_edges.insert({x, y, new_type});

        // We have updated edge (x, y)
        // We need to check the adjacent edges that might be affected by this update
        add_edges_around(x, y, edges_to_check, new_type == EdgeType::DIRECTED);
    }
}


// insert
// i like the general maintain_cpdag function;
// i should keep track of v-structured edges: not for now

std::ostream &operator<<(std::ostream &os, const PDAG &obj) {
    os << "PDAG: directed edges = {";
    for (auto edge: obj.get_directed_edges()) { os << "(" << edge.first << "-->" << edge.second << "), "; }
    os << "}, undirected edges = {";
    for (auto edge: obj.get_undirected_edges()) { os << "(" << edge.first << "---" << edge.second << "), "; }
    os << "}";
    return os;
}