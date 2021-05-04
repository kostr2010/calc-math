#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

#include <iostream>
#include <vector>

#include "dsl.h"

// ====================
// MATRIX

class Matrix {
  public:
    NO_DEFAULT_CTOR(Matrix);
    COPY_SEMANTICS(Matrix);
    DEFAULT_MOVE_SEMANTICS(Matrix);
    DEFAULT_DTOR(Matrix);

    Matrix(const size_t n_row, const size_t n_col);
    Matrix(const std::vector<std::vector<double>> m);

    friend Matrix operator*(const Matrix& l, const Matrix& r);
    friend Matrix operator*(const Matrix& l, const double r);
    friend Matrix operator*(const double l, const Matrix& r);
    Matrix& operator*=(const Matrix& r);
    Matrix& operator*=(const double r);

    friend Matrix operator/(const Matrix& l, const double r);
    Matrix& operator/=(const double r);

    friend Matrix operator+(const Matrix& l, const Matrix& r);
    Matrix& operator+=(const Matrix& r);

    friend Matrix operator-(const Matrix& l, const Matrix& r);
    Matrix& operator-=(const Matrix& r);

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

    void Transponse();

    double Norm1() const;
    double Norm2() const;
    double Norm3() const;

    double Determinant() const;
    Matrix Minor(const size_t row, const size_t col) const;
    std::vector<double> GetEigenvalues(const double& precision = 1e-12) const;
    // FIXME: add IsTriangular for instant eigenvalues
    double GetMaxEigenvalue(const double& precision = 1e-12) const;

    double& at(const size_t row, const size_t col);
    const double& at(const size_t row, const size_t col) const;

    void ExchangeRow(const size_t row1, const size_t row2);
    void ExchangeCol(const size_t col1, const size_t col2);

    std::pair<size_t, size_t> FindBiggestElement(size_t offset = 0) const;
    std::pair<size_t, size_t> FindSmallestElement(size_t offset = 0) const;

    size_t nRow() const;
    size_t nCol() const;

    bool IsEmpty() const;
    bool IsSquare() const;
    bool IsSymmetrical() const;
    // FIXME: add IsTriangular
    bool IsHessenberg(bool* is_upper_hessenberg = nullptr /* out */) const;

    static Matrix E(const size_t size); // identity matrix
    static Matrix N(const size_t size); // zero matrix
    static Matrix G(const size_t size, const size_t i, const size_t j,
                    const double& phi); // Givens rotation matrix
    static Matrix G(const size_t size, const size_t i, const size_t j, const double& sin,
                    const double& cos); // Givens rotation matrix

    void Nullify();
    void RoundUp(const double precision = 1e-12);

  private:
    std::vector<double> GetEigenvaluesRotation(const double& precision) const;
    std::vector<double> GetEigenvaluesQRHessenberg(const double& precision,
                                                   const bool is_upper_hessenberg) const;
    void QRFactorizationHessenberg(const bool is_upper_hessenberg, Matrix H, Matrix* Q /* out */,
                                   Matrix* R /* out */) const;
    // FIXME: add more eigenvalue algorithms

    size_t n_row_;
    size_t n_col_;
    bool transp_;

    std::vector<std::vector<double>> rows_;
};

#endif