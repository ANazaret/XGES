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
    cnpy::NpyArray arr = cnpy::npy_load("arr2.npy");
    RowMajorMatrixXd m = Eigen::Map<RowMajorMatrixXd>(arr.data<double>(), arr.shape[0], arr.shape[1]);

    BICScorer scorer(m, 1);
    //    std::cout << scorer.local_diff_score(4, {0}, 3) << std::endl;

    XGES xges(m, &scorer);

    xges.fit_heuristic();


    return 0;
}