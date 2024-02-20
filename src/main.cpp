#include "BICScorer.h"
#include "PDAG.h"
#include "XGES.h"
#include "utils.h"

#include <filesystem>
#include <fstream>
#include <iostream>

#include "cnpy/cnpy.h"
#include "cxxopts.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace fs = std::filesystem;
using namespace std::chrono;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        RowMajorMatrixXd;

template<typename T>
RowMajorMatrixXd load_npy(const std::string &filename) {
    cnpy::NpyArray arr = cnpy::npy_load(filename);
    return Eigen::Map<RowMajorMatrixXd>(arr.data<T>(), arr.shape[0], arr.shape[1]);
}

int main(int argc, char *argv[]) {
    srand(0);
    std::cout << std::setprecision(16);
    // Set up the logger
    auto logger = spdlog::stdout_color_mt("stdout_logger");
    logger->set_pattern("[%^%l%$] %v");

    // Define the command line options
    cxxopts::Options options("xges", "Run XGES on given data");
    auto option_adder = options.add_options();
    option_adder("input", "Input data numpy file (must be a contiguous C array)",
                 cxxopts::value<std::string>());
    option_adder("output", "Output file (default `xges-graph.txt`)",
                 cxxopts::value<std::string>()->default_value("xges-graph.txt"));
    option_adder("alpha,a", "Alpha parameter for the BIC score (default 2.)",
                 cxxopts::value<double>()->default_value("2."));
    option_adder("stats", "File to save statistics (default `xges-stats.txt`)",
                 cxxopts::value<std::string>()->default_value("xges-stats.txt"));
    option_adder("0,xges0", "Do not perform the extended search of XGES, just XGES-0.",
                 cxxopts::value<bool>()->default_value("false"));
    option_adder("graph_truth,g", "Graph truth file", cxxopts::value<std::string>());
    option_adder("verbose,v", "Level of verbosity (0-3)",
                 cxxopts::value<int>()->default_value("1"));
    auto args = options.parse(argc, argv);

    if (int verbose_level = args["verbose"].as<int>(); verbose_level == 0) {
        logger->set_level(spdlog::level::off);
    } else if (verbose_level == 1) {
        logger->set_level(spdlog::level::info);
    } else if (verbose_level == 2) {
        logger->set_level(spdlog::level::debug);
    } else if (verbose_level == 3) {
        logger->set_level(spdlog::level::trace);
    } else {
        throw std::runtime_error("Invalid verbose level");
    }

    // Parse the command line options
    fs::path data_path = args["input"].as<std::string>();
    fs::path output_path = args["output"].as<std::string>();
    double alpha = args["alpha"].as<double>();

    logger->info("Loading input: {}", data_path);
    RowMajorMatrixXd m;
    if (data_path.extension() == ".npy") {
        m = load_npy<double>(data_path);
    } else {
        throw std::runtime_error("Input file must be a .npy file");
    }
    int n_variables = m.cols();
    int n_samples = m.rows();
    logger->info("Input loaded. Shape: {} x {}", n_samples, n_variables);
    logger->debug("m[0, 0:2] = {} {}", m(0, 0), m(0, 1));

    logger->info("Computing covariance matrix");
    auto start_time = high_resolution_clock::now();
    BICScorer scorer(m, alpha);
    double elapsed_secs = measure_time(start_time);
    logger->info("Covariance computed in {} seconds", elapsed_secs);
    m.resize(0, 0);// free the memory of m

    XGES xges(n_variables, &scorer);

    if (args.count("graph_truth") > 0) {
        auto ground_truth_pdag = std::make_unique<PDAG>(
                PDAG::from_file(args["graph_truth"].as<std::string>()));
        xges.ground_truth_pdag = std::move(ground_truth_pdag);
    }

    bool extended_search = !args["0"].as<bool>();
    logger->info("Fitting XGES with extended search: {}", extended_search);
    start_time = high_resolution_clock::now();
    xges.fit_xges(extended_search);
    elapsed_secs = measure_time(start_time);

    logger->info("XGES search completed in {} seconds", elapsed_secs);
    logger->info("Score: {}", xges.get_score());

    // Save the output
    std::ofstream out_file(output_path);
    out_file << xges.get_pdag().get_adj_string();
    out_file.close();

    // Save the statistics
    std::ofstream stats_file(args["stats"].as<std::string>());
    stats_file << std::setprecision(16);
    stats_file << "time, " << elapsed_secs << std::endl;
    stats_file << "score, " << xges.get_score() << std::endl;
    stats_file << "score check, " << scorer.score_pdag(xges.get_pdag()) << std::endl;
    stats_file << "score_empty, " << xges.get_initial_score() << std::endl;
    stats_file << "score_increase, " << xges.get_score() - xges.get_initial_score()
               << std::endl;
    for (auto &kv: xges.statistics) {
        stats_file << kv.first << ", " << kv.second << std::endl;
    }
    for (auto &kv: xges.get_pdag().statistics) {
        stats_file << kv.first << ", " << kv.second << std::endl;
    }
    for (auto &kv: scorer.statistics) {
        stats_file << kv.first << ", " << kv.second << std::endl;
    }
    return 0;
}


void test_pdag() {
    PDAG pdag(10);
    pdag.add_undirected_edge(0, 1);
    pdag.add_undirected_edge(1, 2);

    PDAG dag_extension_true1(10);
    dag_extension_true1.add_directed_edge(0, 1);
    dag_extension_true1.add_directed_edge(1, 2);
    PDAG dag_extension_true2(10);
    dag_extension_true2.add_directed_edge(2, 1);
    dag_extension_true2.add_directed_edge(1, 0);
    assert(pdag.get_dag_extension() == dag_extension_true1 ||
           pdag.get_dag_extension() == dag_extension_true2);
}