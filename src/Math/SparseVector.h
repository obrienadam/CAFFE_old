#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H

#include <petsc.h>

class SparseMatrix;

class SparseVector
{
private:

    Vec vec_;

    PetscErrorCode errorCode_;

    // MPI related
    int iLower_, iUpper_;

public:

    SparseVector();
    SparseVector(int m);
    ~SparseVector();

    void allocate(int m);
    void allocate(const SparseVector& other);
    void deallocate();
    void setValue(int i, double value);
    void zeroEntries();
    double operator()(int i);

    double l1Norm();
    double l2Norm();
    double infNorm();

    void print();

    friend class SparseMatrix;
    friend void multiply(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &C);
    friend void multiply(const SparseMatrix &A, const SparseVector &x, SparseVector &b);
    friend void multiplyAdd(const SparseMatrix& A, const SparseVector& x1, const SparseVector& x2, SparseVector& b);
    friend void scale(double alpha, SparseVector& vec);
};

#endif
