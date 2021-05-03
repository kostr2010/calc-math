#include <cmath>

#include "lib-inl.h"

std::function<diff_sheme::point_t(diff_sheme::point_t)> K = [](diff_sheme::point_t x) {
    return exp(cos(x));
};
std::function<diff_sheme::point_t(diff_sheme::point_t)> Q = [](diff_sheme::point_t x) {
    return exp(sin(x));
};
std::function<diff_sheme::point_t(diff_sheme::point_t)> F = [](diff_sheme::point_t x) {
    return 1.0;
};

// const diff_sheme::point_t K_FIXED = K(0.5);
// const diff_sheme::point_t Q_FIXED = Q(0.5);
// const diff_sheme::point_t F_FIXED = F(0.5);

constexpr size_t N_NODES = 11 + 1;

int main(int argc, char* argv[]) {
    diff_sheme::DiffEquation eq(K, Q, F, {1, 1}, {0, 0});

    eq.Solve({0, 1}, N_NODES);

    return 0;
}
