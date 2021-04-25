#ifndef __LIB_H_INCLUDED__
#define __LIB_H_INCLUDED__

#include <array>
#include <functional>
#include <iostream>

namespace RungeKutta {

template <size_t SIZE>
struct ButcherTable {
    using Line = std::array<float, SIZE>;

    explicit constexpr ButcherTable() = default;
    /**
     * c_0 | a_00| a_01| ... | a_0n|
     * c_1 | a_10| a_11| ... | a_1n|
     * ... | ... | ... | ... | ... |
     * c_n | a_n0| a_n1| ... | a_nn|
     * ----------------- ... ------
     *  *  | b_0 | b_1 | ... | b_n |
     */

    explicit constexpr ButcherTable(const Line c, const std::array<Line, SIZE> a, const Line b)
        : c_(c), a_(a), b_(b) {
    }

    bool IsExplicit() const {
        for (size_t i = 0; i < SIZE; i++) {
            for (size_t j = i; j < SIZE; j++) {
                if (a_[i][j] != 0.0) {
                    return false;
                }
            }
        }

        return true;
    }

    Line c_;
    std::array<Line, SIZE> a_;
    Line b_;
};

constexpr ButcherTable<4> RK4({0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0},
                              {ButcherTable<4>::Line{0, 0, 0, 0},
                               ButcherTable<4>::Line{1.0 / 2.0, 0, 0, 0},
                               ButcherTable<4>::Line{0, 1.0 / 2.0, 0, 0},
                               ButcherTable<4>::Line{0, 0, 1, 0}},
                              {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0});

template <size_t DIM_VAR>
using Vec = std::array<float, DIM_VAR>;

template <size_t DIM_VAR>
using Func = std::function<Vec<DIM_VAR>(const float x, const Vec<DIM_VAR>& y)>;

template <size_t ORDER, size_t DIM_VAR>
Vec<DIM_VAR> CalculateRK(const ButcherTable<ORDER>& table, const Vec<DIM_VAR>& u_0,
                         const Func<DIM_VAR>& f, const float at, const float step);

}; // namespace RungeKutta

#endif