#include <cassert>
#include <iomanip>
#include <iostream>

#include "linear-system.h"

SLAE::SLAE(const Matrix& A, const Matrix& f) : A_(A), f_(f) {
    assert(A.IsSquare());
    assert(!A.IsEmpty());

    assert(!f.IsEmpty());
    assert(f.nCol() == 1);
    assert(f.nRow() == A.nRow());

    size_ = A.nRow();

    for (size_t i = 0; i < size_; i++) {
        relative_order_.push_back(i);
    }
}

SLAE::SLAE(const SLAE& r) : A_(r.A_), f_(r.f_) {
    assert(r.IsValid());

    size_           = r.size_;
    relative_order_ = r.relative_order_;
}

SLAE& SLAE::operator=(const SLAE& r) {
    assert(IsValid());
    assert(r.IsValid());

    size_           = r.size_;
    relative_order_ = r.relative_order_;

    return *this;
}

std::ostream& operator<<(std::ostream& os, const SLAE& s) {
    os << ".--\n";

    for (size_t i = 0; i < s.size_; i++) {
        os << "| |";

        for (size_t j = 0; j < s.size_; j++) {
            os << '(' << std::setw(10) << std::right << s.A_.at(i, j) << " * x"
               << s.relative_order_.at(j) << ')';

            if (j != s.size_ - 1) {
                os << " + ";
            }
        }

        os << "| = |" << s.f_.at(i, 0) << '|' << std::endl;
    }

    os << "*--" << std::endl;

    return os;
}

SLAE operator*(const SLAE& l, const double r) {
    return SLAE(l.A_ * r, l.f_ * r);
}

SLAE operator*(const double l, const SLAE& r) {
    return SLAE(r.A_ * l, r.f_ * l);
}

SLAE& SLAE::operator*=(const double r) {
    this->A_ *= r;
    this->f_ *= r;

    return *this;
}

SLAE operator/(const SLAE& l, const double r) {
    assert(r);

    return SLAE(l.A_ / r, l.f_ / r);
}

SLAE& SLAE::operator/=(const double r) {
    assert(r);

    this->A_ /= r;
    this->f_ /= r;

    return *this;
}

void SLAE::DivRow(const size_t row, const double divisor) {
    assert(row < size_);
    assert(divisor);

    for (size_t j = 0; j < size_; j++) {
        this->A_.at(row, j) /= divisor;
    }

    this->f_.at(row, 0) /= divisor;
}

void SLAE::MulRow(const size_t row, const double multiplicator) {
    assert(row < size_);

    for (size_t j = 0; j < size_; j++) {
        this->A_.at(row, j) *= multiplicator;
    }

    this->f_.at(row, 0) *= multiplicator;
}

void SLAE::ExchangeRow(const size_t row1, const size_t row2) {
    assert(row1 < size_);
    assert(row2 < size_);

    if (row1 == row2) {
        return;
    }

    this->f_.ExchangeRow(row1, row2);
    this->A_.ExchangeRow(row1, row2);
}

void SLAE::ExchangeCol(const size_t col1, const size_t col2) {
    assert(col1 < size_);
    assert(col2 < size_);

    if (col1 == col2) {
        return;
    }

    this->A_.ExchangeCol(col1, col2);

    const auto cpy                 = this->relative_order_.at(col1);
    this->relative_order_.at(col1) = this->relative_order_.at(col2);
    this->relative_order_.at(col2) = cpy;
}

void SLAE::AddRow(const size_t row1, const size_t row2, const double coeff1 /* = 1 */,
                  const double coeff2 /* = 1 */) {
    assert(row1 < size_);
    assert(row2 < size_);

    for (size_t i = 0; i < size_; i++) {
        A_.at(row1, i) = A_.at(row1, i) * coeff1 + A_.at(row2, i) * coeff2;
    }

    f_.at(row1, 0) = f_.at(row1, 0) * coeff1 + f_.at(row2, 0) * coeff2;
}

std::pair<size_t, size_t> SLAE::FindBiggestElement(size_t offset) const {
    return A_.FindBiggestElement(offset);
}

std::pair<size_t, size_t> SLAE::FindSmallestElement(size_t offset) const {
    return A_.FindSmallestElement(offset);
}

bool SLAE::IsValid() const {
    bool f_is_col    = f_.nCol() == 1;
    bool sizes_match = f_.nRow() == A_.nRow();

    return f_is_col && sizes_match && !A_.IsEmpty() && A_.IsSquare() && !f_.IsEmpty();
}

size_t SLAE::GetSize() const {
    return size_;
}

size_t SLAE::GetVarByCol(size_t col) const {
    assert(col < GetSize());

    return relative_order_[col];
}

void SLAE::RoundUp(const double precision) {
    A_.RoundUp(precision);
    f_.RoundUp(precision);

    return;
}