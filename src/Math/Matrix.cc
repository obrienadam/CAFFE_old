/**
 * @file    Matrix.cc
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
 * Contains the implementations for class Matrix, which uses blas and lapacke
 * as the underlying engine.
 */

#include <cstdlib>
#include <iostream>
#include <algorithm>

extern "C"
{
    #include <cblas.h>
    #include <clapack.h>
}

#include "Matrix.h"
#include "Output.h"

Matrix::Matrix(int m, int n)
{
    allocate(m, n);
}

Matrix::Matrix(double **elements, int m, int n)
    :
      Matrix(m, n)
{
    int i, j;

    for(i = 0; i < m_; ++i)
    {
        for(j = 0; j < n_; ++j)
        {
            operator ()(i, j) = elements[i][j];
        }
    }
}

Matrix::Matrix(double *elements, int m, int n)
    :
      Matrix(m, n)
{
    int k;

    for(k = 0; k < nElements_; ++k)
    {
        elements_[k] = elements[k];
    }
}

Matrix::Matrix(const Matrix &other)
    :
      Matrix(other.m_, other.n_)
{
    elements_ = other.elements_;
}

Matrix::~Matrix()
{

}

Matrix& Matrix::operator=(const Matrix& rhs)
{
    if(this == &rhs)
        return *this;

    reallocate(rhs.m_, rhs.n_);
    elements_ = rhs.elements_;

    return *this;
}

void Matrix::allocate(int m, int n)
{
    m_ = m;
    n_ = n;
    nElements_ = m_*n_;

    elements_.resize(nElements_);
    ipiv_.resize(nElements_);
}

void Matrix::reallocate(int m, int n)
{
    allocate(m, n);
}

void Matrix::addVector3DToRow(const Vector3D &vec, int i, int j)
{
    elements_[n_*i + j] = vec.x;
    elements_[n_*i + j + 1] = vec.y;
    elements_[n_*i + j + 2] = vec.z;
}

double& Matrix::operator()(int i, int j)
{
    if(i > m_ - 1 || i < 0 || j > n_ - 1 || j < 0)
        Output::raiseException("Matrix", "operator()", "tried to access element outside of the matrix.");

    return elements_[n_*i + j];
}

void Matrix::solveLeastSquares(Matrix &b)
{
    throw "\"Matrix::solveLeastSquares\" not currently working because of an apparent bug in clapack_dgels.";
}

void Matrix::solve(Matrix &b)
{
    clapack_dgesv(CblasRowMajor, m_, b.n_, elements_.data(), m_, ipiv_.data(), b.elements_.data(), b.m_);
}

Matrix& Matrix::inverse()
{
    clapack_dgetrf(CblasRowMajor, m_, n_, elements_.data(), n_, ipiv_.data());
    clapack_dgetri(CblasRowMajor, m_, elements_.data(), n_, ipiv_.data());
    return *this;
}

void Matrix::print()
{
    int i, j;

    for(i = 0; i < m_; ++i)
    {
        for(j = 0; j < n_; ++j)
        {
            std::cout << operator()(i, j) << " ";
        } // end for j

        std::cout << std::endl;
    } // end for i
}

Matrix solveLeastSquares(Matrix A, Matrix b)
{
    A.solveLeastSquares(b);
    return b;
}

Matrix solve(Matrix A, Matrix b)
{
    A.solve(b);
    return b;
}

Matrix eye(int m)
{
    Matrix id(m, m);

    for(int i = 0; i < m; ++i)
        id(i, i) = 1.;

    return id;
}

Matrix random(int m, int n, double min, double max)
{
    int i, j;

    Matrix mat(m, n);

    for(i = 0; i < mat.nRows(); ++i)
    {
        for(j = 0; j < mat.nCols(); ++j)
        {
            mat(i, j) = min + double(rand())/double(RAND_MAX)*(max - min);
        } // end for j
    } // end for i

    return mat;
}

Matrix inverse(Matrix matrix)
{
    matrix.inverse();
    return matrix;
}

Matrix operator*(const Matrix& A, const Matrix& B)
{
    Matrix C(A.m_, B.n_);

    multiply(A, B, C);

    return C;
}

void multiply(const Matrix &A, const Matrix &B, Matrix &C)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.m_, B.n_, A.n_, 1., A.elements_.data(), A.n_, B.elements_.data(), B.n_, 1., C.elements_.data(), C.n_);
}

