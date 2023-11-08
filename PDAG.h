//
// Created by Achille Nazaret on 11/3/23.
//

#pragma once

#include "EdgeQueueSet.h"
#include "Insert.h"
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

// Define the PDAGModification enum class.
enum class PDAGModification {
    REMOVE_DIRECTED_EDGE = 1,
    REMOVE_UNDIRECTED_EDGE = 2,
    ADD_DIRECTED_EDGE = 3,
    ADD_UNDIRECTED_EDGE = 4
};

// Declare the PDAG class.
class PDAG {
private:
    int num_nodes;
    std::vector<int> nodes;
    std::vector<std::set<int>> children;
    std::vector<std::set<int>> parents;
    std::vector<std::set<int>> neighbors;
    std::vector<std::set<int>> adjacent;
    std::vector<std::set<int>> adjacent_reachable;

    int number_of_undirected_edges = 0;
    int number_of_directed_edges = 0;

    std::vector<int> node_version;
    int graph_version = 0;
    std::vector<std::tuple<PDAGModification, int, int>> modifications_history;

public:
    explicit PDAG(int num_nodes);

    /**
     * Returns the number of edges in the PDAG.
     * @return The number of edges in the PDAG.
     */
    int get_number_of_edges() const;


    int get_node_version(int node) const;

    std::set<std::pair<int, int>> get_directed_edges() const;

    std::set<std::pair<int, int>> get_undirected_edges() const;

    const std::set<int> &get_parents(int node) const;

    const std::set<int> &get_children(int node) const;

    const std::set<int> &get_neighbors(int node) const;

    const std::set<int> &get_adjacent(int node) const;

    const std::set<int> get_adjacent_reachable(int node) const;

    std::set<int> get_neighbors_adjacent(int node_y, int node_x) const;

    bool is_clique(const std::set<int> &nodes) const;

    bool has_directed_edge(int x, int y) const;

    bool has_undirected_edge(int x, int y) const;

    void remove_directed_edge(int x, int y);

    void remove_undirected_edge(int x, int y);

    void add_directed_edge(int x, int y);

    void add_undirected_edge(int x, int y);

    bool is_insert_valid(const Insert &insert) const;

    bool block_semi_directed_paths(int src, int dst, const std::set<int> &blocked_nodes) const;

    void apply_insert(const Insert &insert);

    void maintain_cpdag(EdgeQueueSet &edges_to_check, std::set<Edge> &changed_edges);

    bool is_oriented_by_meek_rule_1(int x, int y) const;

    bool is_oriented_by_meek_rule_2(int x, int y) const;

    bool is_oriented_by_meek_rule_3(int x, int y) const;

    bool is_part_of_v_structure(int x, int y) const;


    const std::vector<int> &get_nodes() const;

    std::set<int> get_neighbors_not_adjacent(int node_y, int node_x) const;

    friend std::ostream &operator<<(std::ostream &os, const PDAG &obj);
};
