#include "BICScorer.h"
#include "PDAG.h"
#include "XGES.h"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace fs = std::filesystem;

#include "dependencies/cxxopts.hpp"


#include "cnpy.h"
#include "test.h"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd;

int main(int argc, char *argv[]) {
    srand(0);// random seed

    cxxopts::Options options("xges", "Run XGES algorithm on a dataset");
    auto option_adder = options.add_options();
    option_adder("input,i", "Input data file", cxxopts::value<std::string>());
    option_adder("output,o", "Output file (default `xges-graph.txt`)",
                 cxxopts::value<std::string>()->default_value("xges-graph.txt"));
    option_adder("alpha,a", "Alpha parameter", cxxopts::value<double>()->default_value("0.5"));
    option_adder("stats", "File to save statistics (default `xges-stats.txt`)",
                 cxxopts::value<std::string>()->default_value("xges-stats.txt"));
    auto result = options.parse(argc, argv);

    fs::path data_path = result["input"].as<std::string>();
    fs::path output_path = result["output"].as<std::string>();
    double alpha = result["alpha"].as<double>();

    cnpy::NpyArray arr = cnpy::npy_load(data_path);
    RowMajorMatrixXd m = Eigen::Map<RowMajorMatrixXd>(arr.data<double>(), arr.shape[0], arr.shape[1]);

    BICScorer scorer(m, alpha);

    XGES xges(m, &scorer);

    clock_t start = clock();
    xges.fit_heuristic();
    clock_t end = clock();
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << elapsed_secs << std::endl;
    std::cout << "Score: " << xges.get_score() << std::endl;
    std::cout << "Score check: " << scorer.score_pdag(xges.get_pdag()) << std::endl;
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

    std::ofstream stats_file(result["stats"].as<std::string>());
    stats_file << "time, " << elapsed_secs << std::endl;
    stats_file << "score, " << xges.get_score() << std::endl;
    stats_file << "score_increase, " << xges.get_score() - xges.initial_score << std::endl;
    for (auto &kv: xges.get_pdag().statistics) { stats_file << kv.first << ", " << kv.second << std::endl; }
    return 0;
}