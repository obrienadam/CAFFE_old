#ifndef MATRIX_H
#define MATRIX_H

#include "Vector3D.h"

class Matrix
{
private:

    int m_, n_, nElements_;
    double* elements_;
    int* ipiv_;

public:

    //- Constructors, copy constructors and destructors

    Matrix(int m = 0, int n = 0);
    Matrix(const Matrix& other);
    ~Matrix();

    Matrix& operator=(const Matrix& rhs);

    //- Memory management

    void allocate(int m, int n);
    void deallocate();

    //- Access

    double& operator()(int i, int j);
    int nRows(){return m_;}
    int nCols(){return n_;}
    int nElements(){return nElements_;}

    void addVector3DToRow(const Vector3D& vec, int i, int j);

    //- Solution methods (lapacke wrappers)

    void solveLeastSquares(Matrix& b);
    void solve(Matrix& b);

    //- Matrix manipulation

    void transpose();

    //- Debugging and printing

    void print();
};

//- Solution functions. These will not modify the original matrices

Matrix solveLeastSquares(Matrix A, Matrix b);
Matrix solve(Matrix A, Matrix b);
Matrix eye(int m);
Matrix random(int m, int n, double min, double max);

#endif
