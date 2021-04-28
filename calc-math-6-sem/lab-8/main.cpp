#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>

#include "lib-inl.h"

constexpr double TAU = 1e-3;

// equation for VIII.11.5
constexpr float SIGMA = 10;
constexpr float B = 8.0 / 3.0;
constexpr float R = 28;

/**
 * u[0] = x
 * u[1] = y
 * u[2] = z
 */
RungeKutta::Vec<3> Equation1(const float t, const RungeKutta::Vec<3>& u) {
    RungeKutta::Vec<3> res{};

    res[0] = -SIGMA * (u[0] - u[1]);         // f1
    res[1] = -u[0] * u[2] + R * u[0] - u[1]; // f2
    res[2] = u[0] * u[1] - B * u[2];         // f3

    return res;
}

// equation for X.9.2
constexpr float a_1 = 1e3;
constexpr float A_1 = 0.5;
constexpr float a_2 = 1e6;
constexpr float A_2 = 10;
constexpr float OMEGA = 1e-2;

/**
 * u[0] = y1
 * u[1] = y2
 */
RungeKutta::Vec<2> Equation2(const float t, const RungeKutta::Vec<2>& u, const float a,
                             const float A) {
    RungeKutta::Vec<2> res{};

    res[0] = a * (-1 * (u[0] * u[0] * u[0] / 3.0 - u[0]) + u[1]); // f1
    res[1] = -1 * u[0] + A * cos(OMEGA * t);                      // f2

    return res;
}

template <size_t DIM_VAR>
using FuncR = std::function<RungeKutta::Vec<DIM_VAR>(
    const float t, const RungeKutta::Vec<DIM_VAR>& u, const float a, const float A)>;

// jacobian for Equation2
std::array<FuncR<2>, 2> J_Equation2 = {[](const float t, const RungeKutta::Vec<2>& u, const float a,
                                          const float A) { // u[0] = df1 / dy1, u[1] = df1 / dy2
                                           RungeKutta::Vec<2> res{};

                                           res[0] = -1 * a * u[0] * u[0] + a;
                                           res[1] = a;

                                           return res;
                                       },
                                       [](const float t, const RungeKutta::Vec<2>& u, const float a,
                                          const float A) { // u[0] = df2 / dy1, u[1] = df2 / dy2
                                           RungeKutta::Vec<2> res{};

                                           res[0] = -1;
                                           res[1] = 0;

                                           return res;
                                       }};

RungeKutta::Vec<2> CalculateRosenbrock(const FuncR<2>& f, const RungeKutta::Vec<2>& u,
                                       const float t, const float step, const float a,
                                       const float A);

int main() {
    double t = 0;

    RungeKutta::Func<3> func = Equation1;

    std::ofstream file{"res_explicit.csv"};
    if (file.bad()) {
        std::cerr << "couldn't open file\n";
        return 1;
    }
    file << "t,x,y,z\n";

    // calculating VIII.11.5
    RungeKutta::Vec<3> u = {1, 1, 1};
    for (t = 0; t <= 50; t += TAU) {
        u = RungeKutta::CalculateRK<4, 3>(RungeKutta::RK4, u, func, t, TAU);

        file << t << ',' << u[0] << ',' << u[1] << ',' << u[2] << '\n';
    }
    file.close();

    file.open("res_implicit.csv");
    if (file.bad()) {
        std::cerr << "couldn't open file\n";
        return 1;
    }
    file << "t,y1_11,y2_11,y1_12,y2_12,y1_21,y2_21,y1_22,y2_22\n";

    // calculating X.9.2
    RungeKutta::Vec<2> v1 = {2, 0};
    RungeKutta::Vec<2> v2 = {2, 0};
    RungeKutta::Vec<2> v3 = {2, 0};
    RungeKutta::Vec<2> v4 = {2, 0};
    for (t = 0; t <= 200; t += TAU) {
        v1 = CalculateRosenbrock(Equation2, v1, t, TAU, a_1, A_1);
        v2 = CalculateRosenbrock(Equation2, v2, t, TAU, a_1, A_2);
        v3 = CalculateRosenbrock(Equation2, v3, t, TAU, a_2, A_1);
        v4 = CalculateRosenbrock(Equation2, v4, t, TAU, a_2, A_2);

        file << t << ',' << v1[0] << ',' << v1[1] << ',' << v2[0] << ',' << v2[1] << ',' << v3[0]
             << ',' << v3[1] << ',' << v4[0] << ',' << v4[1] << '\n';
    }
    file.close();

    return 0;
}

std::array<RungeKutta::Vec<2>, 2> operator*(const std::array<RungeKutta::Vec<2>, 2>& l,
                                            const std::array<RungeKutta::Vec<2>, 2>& r) {
    std::array<RungeKutta::Vec<2>, 2> res = {};

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            res[i][j] = l[i][0] * r[0][j] + l[i][1] * r[1][j];
        }
    }

    return res;
}

std::array<RungeKutta::Vec<2>, 2> operator*(const double l,
                                            const std::array<RungeKutta::Vec<2>, 2>& r) {
    std::array<RungeKutta::Vec<2>, 2> res = r;

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            res[i][j] *= l;
        }
    }

    return res;
}

std::array<RungeKutta::Vec<2>, 2> operator*(const std::array<RungeKutta::Vec<2>, 2>& l,
                                            const double r) {
    std::array<RungeKutta::Vec<2>, 2> res = l;

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            res[i][j] *= r;
        }
    }

    return res;
}

std::array<RungeKutta::Vec<2>, 2> operator-(const std::array<RungeKutta::Vec<2>, 2>& l,
                                            const std::array<RungeKutta::Vec<2>, 2>& r) {
    std::array<RungeKutta::Vec<2>, 2> res = l;

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            res[i][j] -= r[i][j];
        }
    }

    return res;
}

std::array<RungeKutta::Vec<2>, 2> operator+(const std::array<RungeKutta::Vec<2>, 2>& l,
                                            const std::array<RungeKutta::Vec<2>, 2>& r) {
    std::array<RungeKutta::Vec<2>, 2> res = l;

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            res[i][j] += r[i][j];
        }
    }

    return res;
}

RungeKutta::Vec<2> operator-(const RungeKutta::Vec<2>& l, const RungeKutta::Vec<2>& r) {
    RungeKutta::Vec<2> res = l;

    for (size_t i = 0; i < 2; i++) {
        res[i] -= r[i];
    }

    return res;
}

RungeKutta::Vec<2> operator+(const RungeKutta::Vec<2>& l, const RungeKutta::Vec<2>& r) {
    RungeKutta::Vec<2> res = l;

    for (size_t i = 0; i < 2; i++) {
        res[i] += r[i];
    }

    return res;
}

RungeKutta::Vec<2> operator*(const RungeKutta::Vec<2>& l, const float r) {
    RungeKutta::Vec<2> res = l;

    for (size_t i = 0; i < 2; i++) {
        res[i] *= r;
    }

    return res;
}

RungeKutta::Vec<2> operator*(const float l, const RungeKutta::Vec<2>& r) {
    RungeKutta::Vec<2> res = r;

    for (size_t i = 0; i < 2; i++) {
        res[i] *= l;
    }

    return res;
}

RungeKutta::Vec<2> operator*(const std::array<RungeKutta::Vec<2>, 2>& l,
                             const RungeKutta::Vec<2>& r) {
    RungeKutta::Vec<2> res = {};

    for (size_t i = 0; i < 2; i++) {
        res[i] = l[i][0] * r[0] + l[i][1] * r[1];
    }

    return res;
}

RungeKutta::Vec<2> CalculateRosenbrock(const FuncR<2>& f, const RungeKutta::Vec<2>& u,
                                       const float t, const float step, const float a,
                                       const float A) {
    static const float rosenbrock_a = 1.077;
    static const float rosenbrock_b = -0.372;
    static const float rosenbrock_c = -0.577;

    static const std::array<RungeKutta::Vec<2>, 2> E = {RungeKutta::Vec<2>{1.0, 0.0},
                                                        RungeKutta::Vec<2>{0.0, 1.0}};

    std::array<RungeKutta::Vec<2>, 2> B = {};

    for (size_t i = 0; i < 2; i++) {
        B[i] = J_Equation2[i](t, u, a, A);
    }

    RungeKutta::Vec<2> f_iter = Equation2(t, u + rosenbrock_c * step * Equation2(t, u, a, A), a, A);

    std::array<RungeKutta::Vec<2>, 2> matr =
        E - rosenbrock_a * step * B - rosenbrock_b * step * step * B * B;

    std::array<RungeKutta::Vec<2>, 2> matr_inv = {};
    matr_inv[0][0] = matr[1][1];
    matr_inv[1][1] = matr[0][0];
    matr_inv[0][1] = -matr[0][1];
    matr_inv[1][0] = -matr[1][0];
    matr_inv = matr_inv * (1.0 / (matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1]));

    return u + step * matr_inv * f_iter;
}