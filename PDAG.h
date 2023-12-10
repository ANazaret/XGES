//
// Created by Achille Nazaret on 11/3/23.
//

#pragma once

#include "CircularBuffer.h"
#include "EdgeQueueSet.h"
#include "Operators.h"
#include "set_ops.h"
#include <map>
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
    const int num_variables;
    const int num_interventions;
    std::vector<int> nodes_variables;
    std::vector<int> nodes_interventions;
    std::vector<int> nodes_all;
    std::vector<FlatSet> children;
    std::vector<FlatSet> parents;
    std::vector<FlatSet> neighbors;
    std::vector<FlatSet> adjacent;
    std::vector<FlatSet> adjacent_reachable;

    int number_of_undirected_edges = 0;
    int number_of_directed_edges = 0;

    std::vector<int> node_version;
    int graph_version = 0;
    std::vector<std::tuple<PDAGModification, int, int>> modifications_history;

    std::vector<char> _block_semi_directed_path_visited;
    std::vector<char> _block_semi_directed_path_blocked;
    CircularBuffer<int> _block_semi_directed_path_queue;

public:
    std::map<std::string, double> statistics;

    explicit PDAG(int num_nodes, int num_interventions = 0);

    // read from file
    static PDAG from_file(const std::string &filename);

    /**
     * Returns the number of edges in the PDAG.
     * @return The number of edges in the PDAG.
     */
    int get_number_of_edges() const;


    int get_node_version(int node) const;

    std::vector<std::pair<int, int>> get_directed_edges() const;
    std::vector<std::pair<int, int>> get_directed_edges_interventions() const;
    std::vector<std::pair<int, int>> get_undirected_edges() const;


    const FlatSet &get_parents(int node) const;

    const FlatSet &get_children(int node) const;

    const FlatSet &get_neighbors(int node) const;

    const FlatSet &get_adjacent(int node) const;

    const FlatSet &get_adjacent_reachable(int node) const;

    FlatSet get_neighbors_adjacent(int node_y, int node_x) const;

    bool is_clique(const FlatSet &nodes) const;

    bool has_directed_edge(int x, int y) const;

    bool has_undirected_edge(int x, int y) const;

    void remove_directed_edge(int x, int y);

    void remove_undirected_edge(int x, int y);

    void add_directed_edge(int x, int y);

    void add_undirected_edge(int x, int y);

    bool is_insert_valid(const Insert &insert, bool reverse = false);

    bool is_reverse_valid(const Reverse &reverse);

    bool is_delete_valid(const Delete &delet) const;

    bool block_semi_directed_paths(int src, int dst, const FlatSet &blocked_nodes,
                                   bool ignore_direct_edge = false);

    void apply_insert(const Insert &insert, EdgeModificationsMap &edge_modifications_map);

    void apply_reverse(const Reverse &reverse, EdgeModificationsMap &edge_modifications_map);

    void apply_delete(const Delete &delet, EdgeModificationsMap &edge_modifications_map);

    void maintain_cpdag(EdgeQueueSet &edges_to_check, EdgeModificationsMap &edge_modifications_map);

    bool is_oriented_by_meek_rule_1(int x, int y) const;

    bool is_oriented_by_meek_rule_2(int x, int y) const;

    bool is_oriented_by_meek_rule_3(int x, int y) const;

    bool is_oriented_by_meek_rule_4(int x, int y) const;

    bool is_part_of_v_structure(int x, int y) const;


    const std::vector<int> &get_nodes_variables() const;

    FlatSet get_neighbors_not_adjacent(int node_y, int node_x) const;

    friend std::ostream &operator<<(std::ostream &os, const PDAG &obj);

    void add_edges_around(int x, int y, EdgeQueueSet &edgeQueueSet, bool is_directed = false);

    PDAG get_dag_extension() const;

    std::string get_adj_string() const;

    inline bool node_is_intervention(int node) const { return node >= num_variables; }

    int shd(const PDAG &other, bool allow_directed_in_other = true) const;
};
