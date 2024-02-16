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
private:
    int n_variables;
    ScorerInterface *scorer;

    PDAG pdag;
    double total_score = 0;
    std::shared_ptr<spdlog::logger> _logger;

    void heuristic_turn_delete_insert(std::vector<Insert> &candidate_inserts,
                                      std::vector<Reverse> &candidate_reverses,
                                      std::vector<Delete> &candidate_deletes,
                                      bool use_threshold, bool initialize_inserts = true);
    void block_each_edge_and_research(bool use_threshold);

    void find_delete_to_y_from_x(int y, int x, const FlatSet &parents_y,
                                 std::vector<Delete> &candidate_deletes, double threshold,
                                 bool directed_xy) const;

public:
    const double initial_score = 0;
    std::unique_ptr<PDAG> ground_truth_pdag;

    XGES(int n_variables, ScorerInterface *scorer);
    XGES(const XGES &other);


    void fit_xges(bool use_threshold, bool extended_search);

    double get_score() const;

    const PDAG &get_pdag() const;

    void find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts,
                           int parent_x = -1, bool positive_only = true);

    void find_deletes_to_y(int y, std::vector<Delete> &candidate_deletes,
                           double threshold = 0) const;

    void find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses);

    void find_reverse_from_x(int x, std::vector<Reverse> &candidate_reverses);

    std::map<std::string, double> statistics;

    double deletion_threshold = -1;

    void update_operator_candidates(EdgeModificationsMap &edge_modifications,
                                    std::vector<Insert> &candidate_inserts,
                                    std::vector<Reverse> &candidate_reverses,
                                    std::vector<Delete> &candidate_deletes);
};
