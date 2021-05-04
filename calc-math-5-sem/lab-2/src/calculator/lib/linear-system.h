#ifndef __LINEAR_SYSTEM_H_INCLUDED__
#define __LINEAR_SYSTEM_H_INCLUDED__

#include <vector>

#include "dsl.h"
#include "matrix.h"

// ====================
// LINEAR SYSTEM

class SLAE {
  public:
    NO_DEFAULT_CTOR(SLAE);
    COPY_SEMANTICS(SLAE);
    DEFAULT_MOVE_SEMANTICS(SLAE);
    DEFAULT_DTOR(SLAE);

    SLAE(const Matrix& A, const Matrix& f);

    friend std::ostream& operator<<(std::ostream& os, const SLAE& s);

    friend SLAE operator*(const SLAE& l, const double r);
    friend SLAE operator*(const double l, const SLAE& r);
    SLAE& operator*=(const double r);

    friend SLAE operator/(const SLAE& l, const double r);
    SLAE& operator/=(const double r);

    void DivRow(const size_t row, const double divisor);
    void MulRow(const size_t row, const double multiplicator);

    void ExchangeRow(const size_t row1, const size_t row2);
    void ExchangeCol(const size_t col1, const size_t col2);

    void AddRow(const size_t row1, const size_t row2, const double coeff1 = 1,
                const double coeff2 = 1); // row1 * coeff1 += row2 * coeff2

    std::pair<size_t, size_t> FindBiggestElement(size_t offset) const;
    std::pair<size_t, size_t> FindSmallestElement(size_t offset) const;

    bool IsValid() const;

    size_t GetSize() const;
    size_t GetVarByCol(size_t col) const;

    void RoundUp(const double precision);

    Matrix A_;
    Matrix f_;

  private:
    size_t size_;

    std::vector<size_t> relative_order_;
};

#endif