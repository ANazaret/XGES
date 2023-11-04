//
// Created by Achille Nazaret on 11/3/23.
//
#include "PDAG.h"

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
    const auto &children_set = children.at(x);
    return children_set.find(y) != children_set.end();
}

bool PDAG::has_undirected_edge(int x, int y) const {
    const auto &neighbors_set = neighbors.at(x);
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
