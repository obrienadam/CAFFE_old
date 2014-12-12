#include <cstdlib>
#include <iostream>

// The following code is to work around a bug on lapack-3.4.2 on Fedora 20, and should be removed once the package
// is updated

//#include <complex.h>
//#define lapack_complex_float    float _Complex
//#define lapack_complex_double   double _Complex

#include <lapacke/lapacke.h>

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

void Matrix::addVector3DToRow(const Vector3D &vec, int i, int j)
{

    elements_[n_*i + j] = vec.x;
    elements_[n_*i + j + 1] = vec.y;
    elements_[n_*i + j + 2] = vec.z;

}

double& Matrix::operator()(int i, int j)
{

    if(i > m_ - 1 || i < 0 || j > n_ - 1 || j < 0)
        Output::raiseException("Matrix", "operator()", "tried to access element outside of the array bounds.");

    return elements_[n_*i + j];

}

void Matrix::solveLeastSquares(Matrix &b)
{

    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m_, n_, b.n_, elements_, n_, b.elements_, b.n_);

    // Modify the dimensions of vector b to reflect the number of unknowns. Note that this doesn't release memory,
    // but it shouldn't matter

    b.m_ = n_;

}

void Matrix::solve(Matrix &b)
{

    LAPACKE_dgesv(LAPACK_ROW_MAJOR, m_, b.n_, elements_, n_, ipiv_, b.elements_, b.n_);

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
