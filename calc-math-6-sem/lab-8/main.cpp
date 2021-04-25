#include <fstream>
#include <iostream>

#include "lib-inl.h"

constexpr float SIGMA = 10;
constexpr float B = 8.0 / 3.0;
constexpr float R = 28;

/**
 * u[0] = x
 * u[1] = y
 * u[2] = z
 */
RungeKutta::Vec<3> Equation(const float t, const RungeKutta::Vec<3>& u) {
    RungeKutta::Vec<3> res{};

    res[0] = -SIGMA * (u[0] - u[1]);
    res[1] = -u[0] * u[2] + R * u[0] - u[1];
    res[2] = u[0] * u[1] - B * u[2];

    return res;
}

int main() {
    double t = 0;
    double tau = 1e-3;

    RungeKutta::Func<3> func = Equation;
    std::ofstream file{"res.csv"};

    if (file.bad()) {
        std::cerr << "couldn't open file\n";
        return 1;
    }

    file << "t,x,y,z\n";

    RungeKutta::Vec<3> u = {1, 1, 1};
    for (t = 0; t <= 50; t += tau) {
        u = RungeKutta::CalculateRK<4, 3>(RungeKutta::RK4, u, func, t, tau);

        file << t << ',' << u[0] << ',' << u[1] << ',' << u[2] << '\n';
    }

    file.close();

    return 0;
}