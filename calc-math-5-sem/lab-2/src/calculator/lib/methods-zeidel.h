#ifndef __ZEIDEL_H_INCLUDED__
#define __ZEIDEL_H_INCLUDED__

#include "linear-system.h"
#include "methods.h"

#include <assert.h>
#include <random>

namespace methods {
class Zeidel {
  public:
    static Matrix Solve(SLAE slae, const double precision = 1e-12, Matrix* residual = nullptr) {
        Matrix L(Zeidel::L(slae.A_));
        Matrix D(Zeidel::D(slae.A_));
        Matrix U(Zeidel::U(slae.A_));

        size_t matrix_sz = slae.GetSize();

        Matrix r(matrix_sz, 1);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        std::default_random_engine re;
        for (size_t i = 0; i < matrix_sz; i++) {
            r.at(i, 0) = unif(re);
        }

        Matrix r_cpy(r);

        double deviation = 0.0;
        do {
            r.Nullify();

            for (size_t i = 0; i < matrix_sz; i++) {
                for (size_t j = i + 1; j < matrix_sz; j++) {
                    r.at(i, 0) -= slae.A_.at(i, j) * r_cpy.at(j, 0);
                }

                for (size_t j = 0; j < i; j++) {
                    r.at(i, 0) -= slae.A_.at(i, j) * r.at(j, 0);
                }

                r.at(i, 0) += slae.f_.at(i, 0);
                r.at(i, 0) /= slae.A_.at(i, i);
            }

            deviation = 0;
            for (size_t i = 0; i < matrix_sz; i++) {
                const double diff = r.at(i, 0) - r_cpy.at(i, 0);

                deviation += diff * diff;
            }
            deviation = std::sqrt(deviation);

            r_cpy = r;

            std::cout << "current deviation: " << std::scientific << deviation
                      << ", expected precision: " << precision << "\n";

        } while (deviation > precision);

        *residual = slae.A_ * r - slae.f_;

        return r;
    }

  private:
    static Matrix L(const Matrix& A) {
        assert(A.IsSquare());

        const size_t matrix_sz = A.nCol();

        Matrix L(matrix_sz, matrix_sz);

        for (size_t i = 0; i < matrix_sz; i++) {
            for (size_t j = 0; j < i; j++) {
                L.at(i, j) = A.at(i, j);
            }
        }

        return L;
    }

    static Matrix U(const Matrix& A) {
        assert(A.IsSquare());

        const size_t matrix_sz = A.nCol();

        Matrix U(matrix_sz, matrix_sz);

        for (size_t i = 0; i < matrix_sz; i++) {
            for (size_t j = i + 1; j < matrix_sz; j++) {
                U.at(i, j) = A.at(i, j);
            }
        }

        return U;
    }

    static Matrix D(const Matrix& A) {
        assert(A.IsSquare());

        const size_t matrix_sz = A.nCol();

        Matrix D(matrix_sz, matrix_sz);

        for (size_t i = 0; i < matrix_sz; i++) {
            D.at(i, i) = A.at(i, i);
        }

        return D;
    }
};
}; // namespace methods

#endif