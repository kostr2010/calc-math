#include <cmath>
#include <fstream>
#include <iostream>

#include "lib-inl.h"

diff_sheme::point_t K(const diff_sheme::point_t x) {
    return exp(cos(x));
}

diff_sheme::point_t Q(const diff_sheme::point_t x) {
    return exp(sin(x));
}

diff_sheme::point_t F(const diff_sheme::point_t x) {
    return 1.0;
}

constexpr size_t N_NODES = 11 + 1;

int main(int argc, char* argv[]) {
    diff_sheme::DiffEquation eq(K, Q, F, {1, 1}, {0, 0});

    std::ofstream out("res.txt");

    // auto res = eq.SolveFixedCoeff({0, 1}, N_NODES, 0.5);

    // out << "# solution of differential equation (fixed)\nx,y\n";
    // res.dump(out);

    // res = eq.VerifyFixedCoeff(res, 0.5);

    // out << "\n# residuals of differential equation (fixed)\nx,y\n";
    // res.dump(out);

    auto res = eq.Solve({0, 1}, N_NODES);

    out << "\n# solution of differential equation\nx,y\n";
    res.dump(out);

    res = eq.Verify(res);

    out << "\n# residuals of differential equation\nx,y\n";
    res.dump(out);

    return 0;
}
