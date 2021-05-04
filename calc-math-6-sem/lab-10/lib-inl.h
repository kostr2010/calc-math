#ifndef LIB_INL_H_INCLUDED
#define LIB_INL_H_INCLUDED

#include "lib.h"

#include <cmath>
#include <fstream>

namespace boundary_problem {
template <size_t N_NODES>
std::array<point_t, N_NODES> SolveBoundary(func_t f, func_t f_y, const pair& segment,
                                           const pair& bounds, const double epsilon) {
    point_t t = 1;
    point_t tpp = 20;

    const point_t step = (segment.second - segment.first) / N_NODES;

    size_t iter = 0;

#ifdef PRINT_RESULTS
    std::ofstream out("res_shooting.csv");
    std::vector<std::array<point_t, N_NODES>> table;
#endif

    while (fabs(tpp - t) > epsilon) {
        const auto y_n = SolveKoshi<N_NODES>(f, segment, {bounds.first, t});

        func_grid_t func = [f_y, y_n, segment, step](size_t index, point_t z) {
            return f_y(segment.first + index * step, y_n[index]);
        };

        const auto y_t_N = IterateKoshi<N_NODES>(func, segment, {0, 1});

        auto swap = t;
        t = t - (y_n[N_NODES - 1] - bounds.second) / y_t_N;
        tpp = swap;

        if (iter++ == 1000) {
            std::cerr << "> computational time exceeded\n";
            break;
        }
#ifdef PRINT_RESULTS
        table.push_back(y_n);
#endif
    }

#ifdef PRINT_RESULTS
    out << 'x';
    for (size_t j = 0; j < table.size(); j++) {
        out << ',' << 'y' << j;
    }
    out << '\n';
    for (size_t i = 0; i < N_NODES; i++) {
        out << segment.first + i * step;
        for (size_t j = 0; j < table.size(); j++) {
            out << ',' << table[j][i];
        }
        out << '\n';
    }
#endif

    return SolveKoshi<N_NODES>(f, segment, {bounds.first, tpp});
}

template <size_t N_NODES>
std::array<point_t, N_NODES> SolveKoshi(func_t f, const pair& segment, const pair& bounds) {
    std::array<point_t, N_NODES> y = {};

    const point_t step = (segment.second - segment.first) / N_NODES;

    y[0] = bounds.first;
    y[1] = 0.5 * step * step * f(segment.first, y[0]) + step * bounds.second + y[0];

    for (size_t i = 2; i < N_NODES; i++) {
        y[i] = step * step * f(segment.first + (i - 1) * step, y[i - 1]) + 2 * y[i - 1] - y[i - 2];
    }

    return y;
}

template <size_t N_NODES>
point_t IterateKoshi(func_t f, const pair& segment, const pair& bounds) {
    std::array<point_t, N_NODES> y = {};

    const point_t step = (segment.second - segment.first) / N_NODES;

    y[0] = bounds.first;
    y[1] = 0.5 * step * step * f(segment.first, y[0]) + step * bounds.second + y[0];

    for (size_t i = 2; i < N_NODES; i++) {
        y[i] = step * step * f(segment.first + (i - 1) * step, y[i - 1]) + 2 * y[i - 1] - y[i - 2];
    }

    return y[N_NODES - 1];
}

template <size_t N_NODES>
point_t IterateKoshi(func_grid_t f, const pair& segment, const pair& bounds) {
    std::array<point_t, N_NODES> y = {};

    const point_t step = (segment.second - segment.first) / N_NODES;

    y[0] = bounds.first;
    y[1] = 0.5 * step * step * f(0, y[0]) + step * bounds.second + y[0];

    for (size_t i = 2; i < N_NODES; i++) {
        y[i] = step * step * f(i - 1, y[i - 1]) + 2 * y[i - 1] - y[i - 2];
    }

    return y[N_NODES - 1];
}

}; // namespace boundary_problem

#endif
