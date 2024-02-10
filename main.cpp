#include "BICScorer.h"
#include "PDAG.h"
#include "XGES.h"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace fs = std::filesystem;
using namespace std::chrono;

#include "dependencies/cxxopts.hpp"


#include "cnpy.h"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd;

template<typename T>
RowMajorMatrixXd load_npy(const std::string &filename) {
    cnpy::NpyArray arr = cnpy::npy_load(filename);
    return Eigen::Map<RowMajorMatrixXd>(arr.data<T>(), arr.shape[0], arr.shape[1]);
}

int main(int argc, char *argv[]) {
    srand(0);// random seed

    std::cout << std::setprecision(16);

    cxxopts::Options options("xges", "Run XGES algorithm on a dataset");
    auto option_adder = options.add_options();
    option_adder("input", "Input data numpy file", cxxopts::value<std::string>());
    option_adder("output", "Output file (default `xges-graph.txt`)",
                 cxxopts::value<std::string>()->default_value("xges-graph.txt"));
    option_adder("alpha,a", "Alpha parameter", cxxopts::value<double>()->default_value("0.5"));
    option_adder("stats", "File to save statistics (default `xges-stats.txt`)",
                 cxxopts::value<std::string>()->default_value("xges-stats.txt"));
    option_adder("interventions",
                 "If provided, numpy file with intervention applied to each sample (-1 if none)",
                 cxxopts::value<std::string>()->default_value(""));

    option_adder("t", "threshold delete", cxxopts::value<bool>());
    option_adder("o", "optimization", cxxopts::value<int>()->default_value("0"));

    option_adder("graph_truth,g", "Graph truth file", cxxopts::value<std::string>());


    auto args = options.parse(argc, argv);

    fs::path data_path = args["input"].as<std::string>();
    fs::path output_path = args["output"].as<std::string>();
    double alpha = args["alpha"].as<double>();

    std::cout << "Loading Input: " << data_path << std::endl;
    RowMajorMatrixXd m;
    if (data_path.extension() == ".npy") {
        m = load_npy<double>(data_path);
    } else if (data_path.extension() == ".npz") {
        cnpy::NpyArray arr = cnpy::npz_load(data_path, "gene_expression");
        m = Eigen::Map<RowMajorMatrixXd>(arr.data<double>(), arr.shape[0], arr.shape[1]);
    } else {
        throw std::runtime_error("Unknown file extension");
    }
    std::cout << "Input shape: " << m.rows() << " " << m.cols() << std::endl;
    std::cout << "m[0, 0:2] = " << m(0, 0) << " " << m(0, 1) << std::endl;

    // todo: uniformize how to configure intervention (number, parameters ...)
    // todo: uniformize num_ and n_

    Eigen::VectorXi m_interventions;
    std::vector<FlatSet> interventions_candidate_variables;
    if (args.count("interventions")) {
        fs::path interventions_path = args["interventions"].as<std::string>();
        if (interventions_path.extension() == ".npy") {
            cnpy::NpyArray arr_interventions = cnpy::npy_load(interventions_path);
            m_interventions = Eigen::Map<Eigen::VectorXi>(arr_interventions.data<int32_t>(),
                                                          arr_interventions.shape[0]);
        } else if (interventions_path.extension() == ".npz") {
            cnpy::NpyArray arr_interventions = cnpy::npz_load(interventions_path, "perturbations");
            m_interventions = Eigen::Map<Eigen::VectorXi>(arr_interventions.data<int32_t>(),
                                                          arr_interventions.shape[0]);
        } else {
            throw std::runtime_error("Unknown file extension");
        }

        // For now assume that intervention i targets variable i
        for (int i = 0; i < m.cols(); ++i) { interventions_candidate_variables.push_back({i}); }
    }
    // fix this by using only one constructor with good defaults
    std::cout << "Computing covariance" << std::endl;
    auto start = high_resolution_clock::now();
    BICScorer scorer =
            (args.count("interventions") > 0) ? BICScorer(m, m_interventions, alpha) : BICScorer(m, alpha);
    double elapsed_secs = duration_cast<duration<double>>(high_resolution_clock::now() - start).count();

    std::cout << "Time computing covariance: " << elapsed_secs << std::endl;

    XGES xges = (args.count("interventions") > 0) ? XGES(m, interventions_candidate_variables, &scorer)
                                                  : XGES(m, &scorer);

    // free the memory of m
    m.resize(0, 0);


    PDAG graph_truth = (args.count("graph_truth") > 0)
                               ? PDAG::from_file(args["graph_truth"].as<std::string>())
                               : PDAG(m.cols());
    if (args.count("graph_truth") > 0) {
        std::cout << "Score truth: " << scorer.score_pdag(graph_truth) << std::endl;
        xges.ground_truth_pdag = &graph_truth;
        std::cout << graph_truth << std::endl;
    }
    // todo: handle it not in this hacky way
    if (args.count("t")) {
        xges.deletion_threshold = -1;
    } else {
        xges.deletion_threshold = 1e-10;
    }
    start = high_resolution_clock::now();
    xges.fit_heuristic(args["o"].as<int>());
    elapsed_secs = duration_cast<duration<double>>(high_resolution_clock::now() - start).count();
    std::cout << "Time searching: " << elapsed_secs << std::endl;
    std::cout << "Score: " << xges.get_score() << std::endl;
    std::cout << "Score check: " << scorer.score_pdag(xges.get_pdag()) << std::endl;
    std::cout << "score_increase, " << xges.get_score() - xges.initial_score << std::endl;
    std::cout << std::endl;


    // print all statistics in xges.get_pdag().statistics
    std::cout << "Statistics: " << std::endl;
    for (auto &kv: xges.get_pdag().statistics) { std::cout << kv.first << ": " << kv.second << std::endl; }
    for (auto &kv: xges.statistics) { std::cout << kv.first << ": " << kv.second << std::endl; }

    std::cout << xges.get_pdag();
    // save the pdag in output_path, in the file pdag.txt. create path, adding "/" if needed

    std::ofstream out_file(output_path);
    out_file << xges.get_pdag().get_adj_string();
    out_file.close();

    std::ofstream stats_file(args["stats"].as<std::string>());
    stats_file << std::setprecision(16);
    stats_file << "time, " << elapsed_secs << std::endl;
    stats_file << "score, " << xges.get_score() << std::endl;
    stats_file << "score_increase, " << xges.get_score() - xges.initial_score << std::endl;
    stats_file << "empty_score, " << xges.initial_score << std::endl;
    for (auto &kv: xges.get_pdag().statistics) { stats_file << kv.first << ", " << kv.second << std::endl; }
    return 0;
}