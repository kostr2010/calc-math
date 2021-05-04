#define PRINT_RESULTS

#include "lib-inl.h"

#undef PRINT_RESULTS

#include <cmath>

int main(int argc, char* argv[]) {
    auto f_1 = [](boundary_problem::point_t x, boundary_problem::point_t y) { return exp(y); };
    auto f_2 = [](boundary_problem::point_t x, boundary_problem::point_t y) { return -exp(y); };
    auto f_y1 = [](boundary_problem::point_t x, boundary_problem::point_t y) { return exp(y); };
    auto f_y2 = [](boundary_problem::point_t x, boundary_problem::point_t y) { return -exp(y); };

    const double a = 1;
    const double b = 1;

    constexpr double EPSILON = 1e-5;
    constexpr size_t N_NODES = 100;

    auto res = boundary_problem::SolveBoundary<N_NODES>(f_2, f_y2, {0, 1}, {a, b}, EPSILON);

    return 0;
}
