//
// Created by Achille Nazaret on 11/6/23.
//

#pragma once

#include "PDAG.h"
#include "ScorerInterface.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;

class XGES {
private:
    int n_variables;
    int n_interventions = 0;
    int n_samples;
    ScorerInterface *scorer;

    std::vector<FlatSet> interventions_candidate_variables;
    std::vector<FlatSet> variables_candidate_interventions;

    PDAG pdag;
    double total_score = 0;

    //    void pre_selection

    void heuristic_turn_delete_insert(std::vector<Insert> &candidate_inserts,
                                      std::vector<Reverse> &candidate_reverses,
                                      std::vector<Delete> &candidate_deletes);
    void initialize_fit(MatrixXd &data, ScorerInterface &scorer);

public:
    const double initial_score = 0;
    XGES(const MatrixXd &data, ScorerInterface *scorer);

    XGES(const Eigen::MatrixXd &data, std::vector<FlatSet> interventions_candidate_variables,
         ScorerInterface *scorer);

    void fit_heuristic();

    double get_score() const;

    const PDAG &get_pdag() const;

    void find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts, int parent_x = -1,
                           bool low_parent_only = false, bool positive_only = true);

    void find_deletes_to_y(int y, std::vector<Delete> &candidate_deletes, double threshold = 0);

    void find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses);

    void find_reverse_from_x(int x, std::vector<Reverse> &candidate_reverses);

    std::map<std::string, double> statistics;

    inline bool node_is_intervention(int node) const { return node >= n_variables; }
    double deletion_threshold = -1;// todo: if we enable starting from a pdag, we need to set this differently


    void update_operator_candidates(EdgeModificationsMap &edge_modifications,
                                    std::vector<Insert> &candidate_inserts,
                                    std::vector<Reverse> &candidate_reverses,
                                    std::vector<Delete> &candidate_deletes);
};
