#include <cstdlib>
#include <iostream>

#include "Matrix.h"
#include "Output.h"

Matrix::Matrix(int m, int n)
    :
      elements_(NULL),
      ipiv_(NULL)
{

    allocate(m, n);

}

Matrix::Matrix(const Matrix &other)
    :
      elements_(NULL),
      ipiv_(NULL)
{

    int i;

    allocate(other.m_, other.n_);

    for(i = 0; i < nElements_; ++i)
        elements_[i] = other.elements_[i];

}

Matrix::~Matrix()
{

    deallocate();

}

Matrix& Matrix::operator=(const Matrix& rhs)
{

    int i;

    if(this == &rhs)
        return *this;

    allocate(rhs.m_, rhs.n_);

    for(i = 0; i < nElements_; ++i)
        elements_[i] = rhs.elements_[i];

    return *this;

}

void Matrix::allocate(int m, int n)
{

    deallocate();

    if(m == 0 || n == 0)
        return;

    m_ = m;
    n_ = n;
    nElements_ = m_*n_;
    elements_ = new double[nElements_];
    ipiv_ = new int[m_];

}

void Matrix::deallocate()
{

    if(elements_ == NULL)
        return;

    delete[] elements_;
    elements_ = NULL;
    delete[] ipiv_;
    ipiv_ = NULL;
    nElements_ = m_ = n_ = 0;

}

double& Matrix::operator()(int i, int j)
{

    //if(i > m_ - 1 || i < 0 || j > n_ - 1 || j < 0)
    // Output::raiseException("Matrix", "operator()", "tried to access element outside of the array bounds.");

    return elements_[n_*i + j];

}

void Matrix::solveLeastSquares(Matrix &b)
{

    int info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m_, n_, b.n_, elements_, n_, b.elements_, b.n_);

    //if(info <= 0)
    //Output::raiseException("Matrix", "solveLeastSquares", "an issue occured with Lapacke routine \"LAPACKE_dgels\".");

}

void Matrix::solve(Matrix &b)
{

    int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, m_, b.n_, elements_, n_, ipiv_, b.elements_, b.n_);

    //if(info <= 0)
    //Output::raiseException("Matrix", "solve", "an issue occurred with Lapacke routine \"LAPACKE_dgesv\".");

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
