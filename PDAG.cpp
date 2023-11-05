//
// Created by Achille Nazaret on 11/3/23.
//
#include <queue>
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

int PDAG::get_number_of_edges() const {
    return number_of_directed_edges + number_of_undirected_edges;
}

int PDAG::get_node_version(int node) const {
    return node_version.at(node);
}

// TODO: Find where this is used and potentially change it to std::vector<std::pair<int, int>>.
std::set<std::pair<int, int>> PDAG::get_directed_edges() const {
    std::set<std::pair<int, int>> edges;
    for (const auto &node: nodes) {
        for (int child: children.at(node)) {
            edges.insert({node, child});
        }
    }
    return edges;
}

// TODO: Find where this is used and potentially change it to std::vector<std::pair<int, int>>.
std::set<std::pair<int, int>> PDAG::get_undirected_edges() const {
    std::set<std::pair<int, int>> edges;
    for (const auto &node: nodes) {
        for (int neighbor: neighbors.at(node)) {
            if (node < neighbor) {
                edges.insert({node, neighbor});
            }
        }
    }
    return edges;
}

const std::set<int> &PDAG::get_parents(int node) const {
    return parents[node];
}

const std::set<int> &PDAG::get_children(int node) const {
    return children[node];
}

const std::set<int> &PDAG::get_neighbors(int node) const {
    return neighbors[node];
}

const std::set<int> &PDAG::get_adjacent(int node) const {
    return adjacent[node];
}

const std::set<int> PDAG::get_adjacent_reachable(int node) const {
    return adjacent_reachable.at(node);
}

std::set<int> PDAG::get_neighbors_adjacent(int node_y, int node_x) const {
    std::set<int> result;
    const auto &neighbors_set = get_neighbors(node_y);
    const auto &adjacent_set = get_adjacent(node_x);

    std::set_intersection(neighbors_set.begin(), neighbors_set.end(),
                          adjacent_set.begin(), adjacent_set.end(),
                          std::inserter(result, result.begin()));
    return result;
}

bool PDAG::is_clique(const std::set<int> &nodes_subset) const {
    for (int node1: nodes_subset) {
        const std::set<int> &adjacent_set = get_adjacent(node1);
        for (int node2: nodes_subset) {
            if (node1 != node2 && adjacent_set.find(node2) == adjacent_set.end()) {
                return false;
            }
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
    if (src == dst) {
        return false;
    }
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

        if (reachable.find(dst) != reachable.end()) {
            return false;
        }
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

    // 1. x and y are not adjacent
    auto &adjacent_x = adjacent.at(x);
    if (adjacent_x.find(y) != adjacent_x.end()) {
        return false;
    }

    // 2. T is a subset of Ne(y) \ Ad(x)
    // <=> T is a subset of Ne(y) and T does not intersect Ad(x)
    auto &neighbors_y = neighbors.at(y);
    if (!is_subset(T, neighbors_y) || have_overlap(T, adjacent_x)) {
        return false;
    }

    // 3. [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y) has not changed
    // <=> insert.effective_parents == [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
    auto ne_y_ad_x_T = neighbors_y;
    // intersect effective_parents with Ad(x)
    intersect_in_place(ne_y_ad_x_T, adjacent_x);
    ne_y_ad_x_T.insert(T.begin(), T.end());
    if (!equal_union(insert.effective_parents, ne_y_ad_x_T, parents.at(y))) {
        return false;
    }

    // 4. [Ne(y) ∩ Ad(x)] ∪ T is a clique
    if (!is_clique(ne_y_ad_x_T)) {
        return false;
    }

    // 5. [Ne(y) ∩ Ad(x)] ∪ T block all semi-directed paths from y to x
    if (!block_semi_directed_paths(y, x, ne_y_ad_x_T)) {
        return false;
    }

    return true;
}

