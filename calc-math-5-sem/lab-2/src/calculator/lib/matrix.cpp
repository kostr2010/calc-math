#include "matrix.h"

#include <assert.h>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

Matrix::Matrix(const size_t n_row, const size_t n_col)
    : n_row_(n_row), n_col_(n_col), transp_(false) {
    rows_.resize(n_row);

    for (size_t i = 0; i < n_row; i++) {
        rows_[i].resize(n_col);
    }
}

Matrix::Matrix(const std::vector<std::vector<double>> m)
    : n_row_(m.size()), n_col_(m.at(0).size()), transp_(false) {
    rows_.resize(n_row_);

    for (size_t i = 0; i < n_row_; i++) {
        rows_[i].resize(n_col_);
    }

    for (size_t i = 0; i < n_row_; i++) {
        assert(m[i].size() == n_col_);

        for (size_t j = 0; j < n_col_; j++) {
            rows_[i][j] = m[i][j];
        }
    }
}

Matrix::Matrix(const Matrix& r) : n_row_(r.n_row_), n_col_(r.n_col_), transp_(r.transp_) {
    rows_.resize(n_row_);

    for (size_t i = 0; i < n_row_; i++) {
        rows_[i].resize(n_col_);

        rows_[i] = r.rows_[i];
    }
}

Matrix& Matrix::operator=(const Matrix& r) {
    assert(this->nCol() == r.nCol());
    assert(this->nRow() == r.nRow());

    for (size_t i = 0; i < r.nRow(); i++) {
        for (size_t j = 0; j < r.nCol(); j++) {
            this->at(i, j) = r.at(i, j);
        }
    }

    return *this;
}

Matrix operator*(const Matrix& l, const Matrix& r) {
    assert(l.nCol() == r.nRow());

    size_t size = l.nCol();

    size_t l_n_row = l.nRow();
    size_t r_n_col = r.nCol();

    Matrix res(l_n_row, r_n_col);

    for (size_t i = 0; i < l_n_row; i++) {
        for (size_t j = 0; j < r_n_col; j++) {
            double sum = 0;

            for (size_t k = 0; k < size; k++) {
                sum += l.at(i, k) * r.at(k, j);
            }

            res.at(i, j) = sum;
        }
    }

    return res;
}

Matrix operator*(const Matrix& l, const double r) {
    Matrix res = l;

    for (auto& row : res.rows_) {
        for (auto& elem : row) {
            elem *= r;
        }
    }

    return res;
}

Matrix operator*(const double l, const Matrix& r) {
    Matrix res = r;

    for (auto& row : res.rows_) {
        for (auto& elem : row) {
            elem *= l;
        }
    }

    return res;
}

Matrix& Matrix::operator*=(const Matrix& r) {
    assert(this->nCol() == r.nRow());
    size_t size = this->nCol();

    size_t l_n_row = this->nRow();
    size_t r_n_col = r.nCol();

    Matrix res(l_n_row, r_n_col);

    for (size_t i = 0; i < l_n_row; i++) {
        for (size_t j = 0; j < r_n_col; j++) {
            double sum = 0;

            for (size_t k = 0; k < size; k++) {
                sum += this->at(i, k) * r.at(k, j);
            }

            res.at(i, j) = sum;
        }
    }

    *this = res;

    return *this;
}

Matrix& Matrix::operator*=(const double r) {
    for (auto& row : rows_) {
        for (auto& elem : row) {
            elem *= r;
        }
    }

    return *this;
}

Matrix operator/(const Matrix& l, const double r) {
    Matrix res = l;

    for (auto& row : res.rows_) {
        for (auto& elem : row) {
            elem /= r;
        }
    }

    return res;
}

Matrix& Matrix::operator/=(const double r) {
    for (auto& row : rows_) {
        for (auto& elem : row) {
            elem /= r;
        }
    }

    return *this;
}

Matrix operator+(const Matrix& l, const Matrix& r) {
    assert(l.nCol() == r.nCol());
    assert(l.nRow() == r.nRow());

    Matrix res = l;

    for (size_t i = 0; i < res.nRow(); i++) {
        for (size_t j = 0; j < res.nCol(); j++) {
            res.at(i, j) = l.at(i, j) + r.at(i, j);
        }
    }

    return res;
}

Matrix& Matrix::operator+=(const Matrix& r) {
    assert(this->nCol() == r.nCol());
    assert(this->nRow() == r.nRow());

    for (size_t i = 0; i < this->nRow(); i++) {
        for (size_t j = 0; j < this->nCol(); j++) {
            this->at(i, j) += r.at(i, j);
        }
    }

    return *this;
}

Matrix operator-(const Matrix& l, const Matrix& r) {
    assert(l.nCol() == r.nCol());
    assert(l.nRow() == r.nRow());

    Matrix res = l;

    for (size_t i = 0; i < res.nRow(); i++) {
        for (size_t j = 0; j < res.nCol(); j++) {
            res.at(i, j) = l.at(i, j) - r.at(i, j);
        }
    }

    return res;
}

Matrix& Matrix::operator-=(const Matrix& r) {
    assert(this->nCol() == r.nCol());
    assert(this->nRow() == r.nRow());

    for (size_t i = 0; i < this->nRow(); i++) {
        for (size_t j = 0; j < this->nCol(); j++) {
            this->at(i, j) -= r.at(i, j);
        }
    }

    return *this;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    if (m.IsEmpty()) {
        os << "matrix is empty!";
    } else {
        for (size_t i = 0; i < m.nRow(); i++) {
            for (size_t j = 0; j < m.nCol(); j++) {
                os << '|' << std::setw(10) << std::left << m.at(i, j) << '|';
            }

            os << std::endl;
        }
    }

    os << std::endl;

    return os;
}

void Matrix::Transponse() {
    transp_ = !transp_;

    return;
}

double Matrix::Norm1() const {
    double res = std::numeric_limits<double>::lowest();
    double sum = 0;

    for (size_t i = 0; i < nRow(); i++) {
        for (size_t j = 0; j < nCol(); j++) {
            sum += this->at(i, j);
        }

        if (sum > res) {
            res = sum;
        }

        sum = 0;
    }

    return res;
}

double Matrix::Norm2() const {
    double res = std::numeric_limits<double>::lowest();
    double sum = 0;

    for (size_t i = 0; i < nCol(); i++) {
        for (size_t j = 0; j < nRow(); j++) {
            sum += this->at(i, j);
        }

        if (sum > res) {
            res = sum;
        }

        sum = 0;
    }

    return res;
}

double Matrix::Norm3() const {
    // FIXME:

    return 0;
}

double Matrix::Determinant() const {
    assert(n_col_ == n_row_);

    size_t size = n_col_;

    double sign = 1;
    double d    = 0;

    switch (size) {
    case 0:
        std::cerr << "taking determinant of empty matrix!\n";

        return 0;
    case 1:
        return this->at(0, 0);
    case 2:
        return this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0);
    default:
        for (size_t i = 0; i < size; i++) {
            d += this->at(i, 0) * sign * Minor(i, 0).Determinant();

            sign *= -1;
        }

        return d;
    }

    return 0;
}

Matrix Matrix::Minor(const size_t row, const size_t col) const {
    assert(row < this->nRow());
    assert(col < this->nCol());

    auto elements = rows_;

    elements.erase(elements.begin() + row);
    for (auto& r : elements) {
        r.erase(r.begin() + col);
    }

    Matrix res(elements);

    return res;
}

std::vector<double> Matrix::GetEigenvalues(const double& precision /* = 1e-12*/) const {
    assert(IsSquare());

    bool is_upper_hessenberg = false;

    if (IsSymmetrical()) {
        return GetEigenvaluesRotation(precision);
    }

    if (IsHessenberg(&is_upper_hessenberg)) {
        return GetEigenvaluesQRHessenberg(precision, is_upper_hessenberg);
    }

    return std::vector<double>{};
}

double Matrix::GetMaxEigenvalue(const double& precision /* = 1e-12 */) const {
    assert(IsSquare());

    size_t matrix_sz = this->nCol();

    Matrix r(matrix_sz, 1);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;
    for (size_t i = 0; i < matrix_sz; i++) {
        r.at(i, 0) = unif(re);
    }

    Matrix r_t(r);
    r_t.Transponse();

    double u     = 0.0;
    double u_cpy = 0.0;

    do {
        u_cpy = u;

        u = (r_t * *this * r).at(0, 0) / (r_t * r).at(0, 0);

        Matrix Ar(*this * r);

        r = (Ar) / (Ar).Norm1();

        r_t.Transponse();
        r_t = r;
        r_t.Transponse();

    } while (std::abs(u - u_cpy) > precision);

    return u;
}

double& Matrix::at(const size_t row, const size_t col) {
    assert(row < this->nRow());
    assert(col < this->nCol());

    if (transp_) {
        return rows_[col][row];
    }

    return rows_[row][col];
}

const double& Matrix::at(const size_t row, const size_t col) const {
    assert(row < this->nRow());
    assert(col < this->nCol());

    if (transp_) {
        return rows_[col][row];
    }

    return rows_[row][col];
}

void Matrix::ExchangeRow(const size_t row1, const size_t row2) {
    assert(row1 < nRow());
    assert(row2 < nRow());

    if (row1 == row2) {
        return;
    }

    auto row_cpy = rows_.at(row1);

    rows_.at(row1) = rows_.at(row2);
    rows_.at(row2) = row_cpy;
}

void Matrix::ExchangeCol(const size_t col1, const size_t col2) {
    assert(col1 < nCol());
    assert(col2 < nCol());

    if (col1 == col2) {
        return;
    }

    for (auto& row : rows_) {
        auto elem_cpy = row.at(col1);

        row.at(col1) = row.at(col2);
        row.at(col2) = elem_cpy;
    }
}

size_t Matrix::nRow() const {
    return (transp_) ? (n_col_) : (n_row_);
}

size_t Matrix::nCol() const {
    return (transp_) ? (n_row_) : (n_col_);
}

std::pair<size_t, size_t> Matrix::FindBiggestElement(size_t offset /* = 0 */) const {
    assert(offset < nCol());
    assert(offset < nRow());

    double min = std::numeric_limits<double>::min();
    std::pair<size_t, size_t> res{};

    for (size_t col = offset; col < nCol(); col++) {
        for (size_t row = offset; row < nRow(); row++) {
            double cur = at(row, col);

            if (std::abs(cur) > min) {
                min = std::abs(cur);

                res = {row, col};
            }
        }
    }

    return res;
}

std::pair<size_t, size_t> Matrix::FindSmallestElement(size_t offset /* = 0 */) const {
    assert(offset < nCol());
    assert(offset < nRow());

    double max = std::numeric_limits<double>::max();
    std::pair<size_t, size_t> res{};

    for (size_t col = offset; col < nCol(); col++) {
        for (size_t row = offset; row < nRow(); row++) {
            double cur = at(row, col);

            if (std::abs(cur) < max) {
                max = std::abs(cur);

                res = {row, col};
            }
        }
    }

    return res;
}

bool Matrix::IsEmpty() const {
    return rows_.size() == 0;
}

bool Matrix::IsSquare() const {
    return n_col_ == n_row_;
}

bool Matrix::IsSymmetrical() const {
    if (!IsSquare()) {
        return false;
    } else if (IsEmpty()) {
        return true;
    }

    for (size_t row = 0; row < nRow(); row++) {
        for (size_t col = 0; col < row; col++) {
            if (this->at(row, col) != this->at(col, row)) {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::IsHessenberg(bool* is_upper_hessenberg /* = nullptr */ /* out */) const {
    if (!IsSquare()) {
        return false;
    }

    bool is_upper = true;
    bool is_lower = true;
    bool is_hess  = true;

    for (size_t i = 1; i < nCol(); i++) {
        for (size_t j = 0; j < i - 1; j++) {
            if (this->at(i, j) != 0.0) {
                is_upper = false;
            }

            if (this->at(j, i) != 0.0) {
                is_lower = false;
            }
        }

        is_hess = is_lower || is_upper;

        if (!is_hess) {
            break;
        }
    }

    if (is_upper_hessenberg != nullptr) {
        *is_upper_hessenberg = is_upper;
    }

    return is_hess;
}

Matrix Matrix::E(const size_t size) {
    Matrix res(size, size);

    for (size_t i = 0; i < size; i++) {
        res.at(i, i) = 1;
    }

    return res;
}

Matrix Matrix::N(const size_t size) {
    return Matrix(size, size);
}

Matrix Matrix::G(const size_t matrix_sz, const size_t i, const size_t j, const double& sin,
                 const double& cos) {
    Matrix G = Matrix::E(matrix_sz);

    G.at(i, j) = -1 * sin;
    G.at(j, i) = sin;
    G.at(i, i) = cos;
    G.at(j, j) = cos;

    return G;
}

Matrix Matrix::G(const size_t matrix_sz, const size_t i, const size_t j, const double& phi) {
    return G(matrix_sz, i, j, std::cos(phi), std::sin(phi));
}

void Matrix::Nullify() {
    for (size_t col = 0; col < nCol(); col++) {
        for (size_t row = 0; row < nRow(); row++) {
            at(row, col) = 0;
        }
    }

    return;
}

void Matrix::RoundUp(const double precision /* = 1e-12 */) {
    for (size_t i = 0; i < nRow(); i++) {
        for (size_t j = 0; j < nCol(); j++) {
            if (std::abs(at(i, j)) < precision) {
                at(i, j) = 0;
            }
        }
    }

    return;
}

std::vector<double> Matrix::GetEigenvaluesRotation(const double& precision) const {
    assert(this->IsSymmetrical());

    Matrix A(*this);

    size_t matrix_sz = A.nCol();

    size_t steps      = 0;
    size_t steps_teor = matrix_sz * (matrix_sz - 1) / 2;

    double deviation = 0;
    do {
        // update rotation matrix
        const auto max_elem = A.FindBiggestElement();

        const auto a_ij = A.at(max_elem.first, max_elem.second);
        const auto a_ii = A.at(max_elem.first, max_elem.first);
        const auto a_jj = A.at(max_elem.second, max_elem.second);

        const auto rot_angle = 0.5 * std::atan(2 * a_ij / (a_ii - a_jj));

        Matrix rot(G(matrix_sz, max_elem.first, max_elem.second, rot_angle));

        // rotate
        A = rot * A;
        rot.Transponse();
        A = A * rot;

        // update deviation
        deviation = 0;
        for (size_t row = 0; row < matrix_sz; row++) {
            for (size_t col = 0; col < row; col++) {
                const auto elem = A.at(row, col);

                deviation += elem * elem;
            }
        }

        deviation = std::sqrt(2 * deviation);

        // update steps
        ++steps;
    } while ((deviation > precision) && (steps > matrix_sz * steps_teor));

    std::vector<double> res{};

    for (size_t i = 0; i < matrix_sz; i++) {
        res.push_back(A.at(i, i));
    }

    return res;
}

std::vector<double> Matrix::GetEigenvaluesQRHessenberg(const double& precision,
                                                       const bool is_upper_hessenberg) const {
    assert(IsHessenberg());
    assert(!IsEmpty());

    size_t matrix_sz = nCol();

    Matrix Q(Matrix::N(matrix_sz));
    Matrix R(Matrix::N(matrix_sz));
    Matrix H(*this);

    double deviation     = std::numeric_limits<double>::max();
    double deviation_cpy = std::numeric_limits<double>::max();
    do {
        QRFactorizationHessenberg(is_upper_hessenberg, H, &Q, &R);

        H = R * Q;

        deviation_cpy = deviation;
        deviation     = 0;
        for (size_t i = 0; i < matrix_sz; i++) {
            if (is_upper_hessenberg) {
                for (size_t j = 0; j < i; j++) {
                    const auto elem = H.at(i, j);

                    deviation += elem * elem;
                }
            } else {
                for (size_t j = i + 1; j < matrix_sz; j++) {
                    const auto elem = H.at(i, j);

                    deviation += elem * elem;
                }
            }
        }

        deviation = std::sqrt(deviation);
        // std::cout << "current deviation: " << std::scientific << deviation
        //           << ", expected precision: " << precision << "\n";
    } while ((deviation > precision) && (deviation < deviation_cpy));

    std::vector<double> res{};

    for (size_t i = 0; i < matrix_sz; i++) {
        res.push_back(H.at(i, i));
    }

    return res;
}

void Matrix::QRFactorizationHessenberg(const bool is_upper_hessenberg, Matrix H,
                                       Matrix* Q /* out */, Matrix* R /* out */) const {
    size_t matrix_sz = H.nCol();

    *Q = Matrix::E(matrix_sz);

    for (size_t k = 0; k < matrix_sz - 1; k++) {
        double i = (is_upper_hessenberg) ? (k) : (matrix_sz - 1 - k);
        double j = (is_upper_hessenberg) ? (k + 1) : (matrix_sz - 1 - 1 - k);

        double H_ii = H.at(i, i);
        double H_ji = H.at(j, i);

        const double r_i = std::sqrt(H_ii * H_ii + H_ji * H_ji);
        const double cos = H_ii / r_i;
        const double sin = H_ji / r_i;

        Matrix rot(G(matrix_sz, i, j, sin, cos));

        *Q *= rot;

        rot.Transponse();

        H = rot * H;
    }

    *R = H;

    return;
}