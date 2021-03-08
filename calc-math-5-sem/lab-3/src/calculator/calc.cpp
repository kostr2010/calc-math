#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>

int main(int argc, char* argv[], char* envp[]) {
    const auto func = [](const double& x) { return log(x) / x; };

    // calculated by hand
    const double x_max = exp(1);
    const double y_max = func(x_max);

    // f_1(x) = 2 * ln(x) / y_max
    // f_2(x) = e^(0.5 * y_max * x)
    const auto f_2 = [&y_max](const double& x) { return 2.0 * log(x) / y_max; };
    const auto f_1 = [&y_max](const double& x) { return exp(0.5 * y_max * x); };

    // |f_1'(x)| must be < 1 for method to converge -> x (5.435; inf) => right solution
    // |f_2'(x)| must be < 1 for method to converge -> x (-inf; 9.2) => left solution

    const double epsilon = 1e-3;
    double x_k           = x_max;
    double x_kk          = 0;
    double x_k_cpy       = x_k;

    while (abs(x_kk - x_k_cpy) > epsilon / 2.0) {
        x_kk    = f_1(x_k);
        x_k_cpy = x_k;
        x_k     = x_kk;
    }

    double res_1 = x_kk;

    x_k  = x_max;
    x_kk = 0;

    while (abs(x_kk - x_k_cpy) > epsilon / 2.0) {
        x_kk    = f_2(x_k);
        x_k_cpy = x_k;
        x_k     = x_kk;
    }

    double res_2 = x_kk;

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

    ofs << "first solution: " << res_1 << "\n";
    ofs << "second solution: " << res_2 << "\n";
    ofs << "width: " << abs(res_2 - res_1) << "\n";

    std::cout << "done!\n";

    return 0;
}