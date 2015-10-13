/**
 * @file    Matrix.h
 * @author  Adam O'Brien <obrienadam89@gmail.com>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * Contains the interface for class Matrix, which is used for storing,
 * solving and performing other linear algebra operations.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <initializer_list>

#include "Vector3D.h"

class Matrix
{
public:

    //- Constructors, copy constructors and destructors

    Matrix(int m = 0, int n = 0);
    Matrix(int m, int n, std::initializer_list<double> list);
    Matrix(const Matrix& other);
    ~Matrix();

    Matrix& operator=(const Matrix& rhs);

    void allocate(int m, int n);
    void reallocate(int m, int n);
    void zero();

    double& operator()(int i, int j);
    double* data(){ return elements_.data(); }
    const double* data() const { return elements_.data(); }
    int nRows() const {return m_;}
    int nCols() const {return n_;}
    int nElements(){return elements_.size();}

    /** Add a 3D vector to a row, with the first element at i, j.
     */
    void addVector3DToRow(const Vector3D& vec, int i, int j);

    /** Solve a least squares problem.
     * @param b a matrix of rhs vectors
     */
    void solveLeastSquares(Matrix& b);

    /** Solve the matrix using an LU factorization.
    * @param b a matrix of rhs vectors
    */
    void solve(Matrix& b);

    /**
     * @brief Compute the inverse using an LU factorization.
     * @return A reference to the inverted matrix.
     */
    Matrix& inverse();

    /**
     * @brief Compute the transpose of the matrix.
     * @return A reference to the transposed matrix.
     */
    Matrix& transpose();

    Matrix& operator+=(Matrix& rhs);
    Matrix& operator-=(Matrix& rhs);
    Matrix& operator*=(double alpha);
    Matrix& operator/=(double alpha);

    /** Print the matrix to the console.
     */
    void print();

protected:

    static std::vector<double> buffer_;

    int m_, n_, nElements_;
    std::vector<double> elements_;
    std::vector<int> ipiv_;

    friend void multiply(const Matrix &A, const Matrix &B, Matrix &C);
    friend Matrix operator*(const Matrix& A, const Matrix& B);
};

//- Solution functions. These will not modify the original matrices

Matrix solveLeastSquares(Matrix A, Matrix b);
Matrix solve(Matrix A, Matrix b);
Matrix eye(int m);
Matrix random(int m, int n, double min, double max);
Matrix inverse(Matrix matrix);
Matrix transpose(Matrix matrix);

Matrix operator*(const Matrix& A, const Matrix& B);
void multiply(const Matrix &A, const Matrix &B, Matrix &C);

#endif
