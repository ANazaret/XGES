//
// Created by Achille Nazaret on 11/3/23.
//

#pragma once

#include "CircularBuffer.h"
#include "EdgeQueueSet.h"
#include "Operators.h"
#include "set_ops.h"
#include <map>
#include <unordered_set>
#include <vector>

// shortcut type for storing paths and pairs of nodes
typedef std::unordered_map<std::pair<int, int>, std::unordered_set<std::pair<int, int>>>
        UnblockedPathsMap;


class PDAG {
    int num_variables;
    std::vector<int> nodes_variables;
    std::vector<int> nodes_all;
    std::vector<FlatSet> children;
    std::vector<FlatSet> parents;
    std::vector<FlatSet> neighbors;
    std::vector<FlatSet> adjacent;
    std::vector<FlatSet> adjacent_reachable;

    int number_of_undirected_edges = 0;
    int number_of_directed_edges = 0;

    std::vector<char> _block_semi_directed_path_visited;
    std::vector<char> _block_semi_directed_path_blocked;
    std::vector<int> _block_semi_directed_path_parent;
    CircularBuffer<int> _block_semi_directed_path_queue;

    std::vector<FlatSet> forbidden_insert_parents;
    std::vector<FlatSet> forbidden_insert_children;

public:
    std::map<std::string, double> statistics;

    explicit PDAG(int num_nodes);
    static PDAG from_file(const std::string &filename);

    bool operator==(const PDAG &other) const;

    int get_number_of_edges() const;
    std::vector<std::pair<int, int>> get_directed_edges() const;
    std::vector<std::pair<int, int>> get_undirected_edges() const;
    const std::vector<int> &get_nodes_variables() const;

    const FlatSet &get_parents(int node) const;
    const FlatSet &get_children(int node) const;
    const FlatSet &get_neighbors(int node) const;
    const FlatSet &get_adjacent(int node) const;
    const FlatSet &get_adjacent_reachable(int node) const;

    FlatSet get_neighbors_adjacent(int node_y, int node_x) const;
    FlatSet get_neighbors_not_adjacent(int node_y, int node_x) const;

    bool has_directed_edge(int x, int y) const;
    bool has_undirected_edge(int x, int y) const;
    void remove_directed_edge(int x, int y);
    void remove_undirected_edge(int x, int y);
    void add_directed_edge(int x, int y);
    void add_undirected_edge(int x, int y);
    void apply_edge_modification(EdgeModification edge_modification, bool undo = false);

    bool is_clique(const FlatSet &nodes) const;
    bool is_blocking_semi_directed_paths(const int y, const int x,
                                         const FlatSet &blocked_nodes,
                                         UnblockedPathsMap &unblocked_paths_map,
                                         const bool ignore_direct_edge = false);

    bool is_insert_valid(const Insert &insert, UnblockedPathsMap &unblocked_paths_map,
                         bool reverse = false);
    bool is_reverse_valid(const Reverse &reverse, UnblockedPathsMap &unblocked_paths_map);
    bool is_delete_valid(const Delete &delet) const;

    void apply_insert(const Insert &insert, EdgeModificationsMap &edge_modifications_map);
    void apply_reverse(const Reverse &reverse,
                       EdgeModificationsMap &edge_modifications_map);
    void apply_delete(const Delete &delet, EdgeModificationsMap &edge_modifications_map);

    bool is_oriented_by_meek_rule_1(int x, int y) const;
    bool is_oriented_by_meek_rule_2(int x, int y) const;
    bool is_oriented_by_meek_rule_3(int x, int y) const;
    bool is_oriented_by_meek_rule_4(int x, int y) const;
    bool is_part_of_v_structure(int x, int y) const;
    void complete_cpdag(EdgeQueueSet &edges_to_check,
                        EdgeModificationsMap &edge_modifications_map);
    bool check_is_cpdag();
    void add_adjacent_edges(int x, int y, EdgeQueueSet &edgeQueueSet);

    bool is_insert_forbidden(const int x, const int y) const {
        return forbidden_insert_parents[y].find(x) != forbidden_insert_parents[y].end();
    }
    void add_forbidden_insert(int parent, int y) {
        forbidden_insert_parents[y].insert(parent);
        forbidden_insert_children[parent].insert(y);
    }

    int shd(const PDAG &other, bool allow_directed_in_other = true) const;
    PDAG get_dag_extension() const;

    std::string get_adj_string() const;
    friend std::ostream &operator<<(std::ostream &os, const PDAG &obj);
};
