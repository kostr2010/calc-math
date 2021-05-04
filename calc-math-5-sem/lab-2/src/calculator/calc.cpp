#include "lib/linear-system.h"
#include "lib/methods-gauss.h"
#include "lib/methods-zeidel.h"
#include "lib/methods.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

const size_t SIZE = 10;

int main(int argc, char* argv[], char* envp[]) {
    const double a = 1;

    Matrix A(SIZE, SIZE);
    Matrix f(SIZE, 1);

    for (size_t i = 0; i < SIZE - 1; i++) {
        double b_ij = 1.0 / (i + 1);
        double f_i  = i + 1;

        for (size_t j = 0; j <= i + 1; j++) {
            A.at(i, j) = (i == j) ? (a) : (b_ij);
        }

        f.at(i, 0) = f_i;
    }

    for (size_t j = 0; j < SIZE; j++) {
        A.at(SIZE - 1, j) = (SIZE - 1 == j) ? (a) : (1.0 / SIZE);
    }

    f.at(SIZE - 1, 0) = SIZE;

    SLAE s{A, f};

    std::cout << "done!\n> calculating eigenvalues (full spectrum research)...\n";

    const auto eigenvalues = A.GetEigenvalues(1e-12);

    std::cout << "done!\n> calculating only min / max eigenvalues...\n";

    const auto eigenvalue_max = A.GetMaxEigenvalue(1e-12);
    auto AA                   = methods::Gauss::Inverse(A);
    const auto eigenvalue_min = 1.0 / AA.GetMaxEigenvalue(1e-12);

    std::cout << "done!\n> solving with Gauss method...\n";

    Matrix residual_gauss(s.GetSize(), 1);
    Matrix res_gauss = methods::Gauss::Solve(s, 1e-12, &residual_gauss);

    std::cout << "done!\n> solving with Zeidel method...\n";

    Matrix residual_zeidel(s.GetSize(), 1);
    Matrix res_zeidel = methods::Zeidel::Solve(s, 1e-12, &residual_zeidel);

    std::cout << "done!\n> emitting answers to output/calculator/result.txt...\n";

    std::string path     = argv[0];
    std::string path_rel = "../output/calculator/";
    std::string path_abs = "";

    size_t dir_switch = path.rfind('/');

    if (dir_switch != std::string::npos) {
        path_abs = path.substr(0, dir_switch + 1) + path_rel;
    } else {
        path_abs = "./" + path_rel;
    }

    std::ofstream ofs(path_abs + "res.txt");

    ofs << "// ====================\n// COMMON\n\n";

    ofs << "eigenvalues for A:\n>";
    for (const auto& elem : eigenvalues) {
        ofs << " " << elem;
    }

    ofs << "\nlambda_min: " << eigenvalue_min << "\n";
    ofs << "lambda_max: " << eigenvalue_max << "\n";

    // lambda_min: 9.34566
    // lambda_max: 11

    ofs << "\n// ====================\n// GAUSS\n";

    ofs << "solution: \n" << res_gauss << "\n";
    ofs << "residual: \n" << residual_gauss << "\n";

    ofs << "\n// ====================\n// ZEIDEL\n";

    ofs << "solution: \n" << res_zeidel << "\n";
    ofs << "residual: \n" << residual_zeidel << "\n";

    std::cout << "done!\n";

    return 0;
}