#ifndef __LIB_INL_H_INCLUDED__
#define __LIB_INL_H_INCLUDED__

#include "lib.h"

#include <iostream>

namespace RungeKutta {

template <size_t ORDER, size_t DIM_VAR>
Vec<DIM_VAR> CalculateRK(const ButcherTable<ORDER>& table, const Vec<DIM_VAR>& u_0,
                         const Func<DIM_VAR>& f, const float at, const float step) {
    auto u = u_0;

    if (table.IsExplicit()) {
        std::array<Vec<DIM_VAR>, ORDER> table_f_i{};

        for (size_t stage = 0; stage < ORDER; stage++) {
            for (size_t i = 0; i < ORDER; i++) {
                for (size_t k = 0; k < DIM_VAR; k++) {
                    u[k] += step * table.b_[i] * table_f_i[i][k];
                }
            }

            auto f_i = f(at + table.c_[stage], u);
            for (size_t i = 0; i < DIM_VAR; i++) {
                table_f_i[stage][i] = f_i[i];
            }
        }
    } else {
        // won't do because it's hard to do general Newton's method
        std::cout << "to be implemented\n";
    }

    return u;
}
}; // namespace RungeKutta

#endif