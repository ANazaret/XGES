//
// Created by Achille Nazaret on 11/6/23.
//

#pragma once

#include "PDAG.h"
#include "ScorerInterface.h"

#define EIGEN_USE_BLAS
#include "Eigen/Dense"
#include "spdlog/logger.h"

using Eigen::MatrixXd;

class XGES {
public:
    XGES(int n_variables, ScorerInterface *scorer);
    XGES(const XGES &other);

    void fit_xges(bool extended_search);
    double get_score() const;
    double get_initial_score() const;
    const PDAG &get_pdag() const;

    std::unique_ptr<PDAG> ground_truth_pdag;
    std::map<std::string, double> statistics;

private:
    int n_variables;
    ScorerInterface *scorer;

    PDAG pdag;
    const double initial_score = 0;
    double total_score = 0;
    std::shared_ptr<spdlog::logger> _logger;

    void heuristic_xges0(std::vector<Insert> &candidate_inserts,
                         std::vector<Reverse> &candidate_reverses,
                         std::vector<Delete> &candidate_deletes,
                         UnblockedPathsMap &unblocked_paths_map,
                         bool initialize_inserts = true);
    void update_operator_candidates_naive(std::vector<Insert> &candidate_inserts,
                                          std::vector<Reverse> &candidate_reverses,
                                          std::vector<Delete> &candidate_deletes);
    void update_operator_candidates_efficient(EdgeModificationsMap &edge_modifications,
                                              std::vector<Insert> &candidate_inserts,
                                              std::vector<Reverse> &candidate_reverses,
                                              std::vector<Delete> &candidate_deletes,
                                              UnblockedPathsMap &unblocked_paths_map);
    void block_each_edge_and_research(UnblockedPathsMap &unblocked_paths_map);

    // todo: make a separate find_inserts_to_y_from_x
    void find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts,
                           int parent_x = -1, bool positive_only = true);

    void find_delete_to_y_from_x(int y, int x, std::vector<Delete> &candidate_deletes,
                                 bool positive_only = true) const;
    void find_deletes_to_y(int y, std::vector<Delete> &candidate_deletes,
                           bool positive_only = true) const;

    void find_reverse_to_y_from_x(int y, int x, std::vector<Reverse> &candidate_reverses);
    void find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses);
    void find_reverse_from_x(int x, std::vector<Reverse> &candidate_reverses);
};
