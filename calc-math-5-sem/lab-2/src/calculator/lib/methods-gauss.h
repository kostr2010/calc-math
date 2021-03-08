#ifndef __GAUSS_H_INCLUDED__
#define __GAUSS_H_INCLUDED__

#include "linear-system.h"
#include "methods.h"

#include <assert.h>

namespace methods {
class Gauss {
  public:
    static Matrix Solve(SLAE slae, const double precision = 1e-12, Matrix* residual = nullptr) {
        Beautify(&slae);
        MakeUpperTriangle(&slae);

        Matrix res(slae.GetSize(), 1);
        res.Nullify();

        SolvePartial(slae, &res, 0);

        if (residual != nullptr) {
            *residual = slae.A_ * res - slae.f_;
        }

        return res;
    }

    static Matrix Inverse(const Matrix& matrix) {
        assert(matrix.IsSquare());

        size_t sz = matrix.nCol();

        Matrix result(sz, sz);
        result.Nullify();

        for (size_t i = 0; i < sz; i++) {
            std::vector<double> part_res(sz, 0);

            part_res[i] = 1;

            Matrix col({part_res});

            col.Transponse();

            SLAE s(matrix, col);

            auto res = Solve(s, 1e-12);
            for (size_t j = 0; j < sz; j++) {
                result.at(j, i) = res.at(j, 0);
            }
        }

        return result;
    }

  private:
    static void SolvePartial(const SLAE& slae, Matrix* res /* out */, size_t offset) {
        size_t var_idx = slae.GetVarByCol(offset);
        assert(slae.A_.at(offset, offset) != 0);

        if (offset + 1 == slae.GetSize()) {
            res->at(var_idx, 0) = slae.f_.at(offset, 0) / slae.A_.at(offset, offset);

            return;
        }

        SolvePartial(slae, res, offset + 1);

        for (size_t i = slae.GetSize() - 1; i > offset; i--) {
            res->at(var_idx, 0) -= res->at(i, 0) * slae.A_.at(offset, i);
        }

        res->at(var_idx, 0) += slae.f_.at(offset, 0);
        res->at(var_idx, 0) /= slae.A_.at(offset, offset);

        return;
    }
};
}; // namespace methods

#endif