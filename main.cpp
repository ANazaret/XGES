#include "BICScorer.h"
#include "PDAG.h"
#include "XGES.h"
#include <iostream>

#include "cnpy.h"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd;

int main() {
    // Seed random number generator of Eigen
    srand(0);
    // read MatrixXd from file
    //save it to file
    cnpy::NpyArray arr = cnpy::npy_load("arr5.npy");
    RowMajorMatrixXd m = Eigen::Map<RowMajorMatrixXd>(arr.data<double>(), arr.shape[0], arr.shape[1]);

    BICScorer scorer(m, 1);

    XGES xges(m, &scorer);

    clock_t start = clock();
    xges.fit_heuristic();
    clock_t end = clock();
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << elapsed_secs << std::endl;
    std::cout << "Score: " << xges.get_score() << std::endl;
    std::cout << std::endl;

    // print all statistics in xges.get_pdag().statistics
    std::cout << "Statistics: " << std::endl;
    for (auto &kv: xges.get_pdag().statistics) { std::cout << kv.first << ": " << kv.second << std::endl; }
    for (auto &kv: xges.statistics) { std::cout << kv.first << ": " << kv.second << std::endl; }

    return 0;
}