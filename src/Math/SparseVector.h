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
    void setValue(int i, double value);
    double operator()(int i);

    void print();

    friend class SparseMatrix;
    friend void multiply(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &C);
    friend void multiply(const SparseMatrix &A, const SparseVector &x, SparseVector &b);
    friend void multiplyAdd(const SparseMatrix& A, const SparseVector& x1, const SparseVector& x2, SparseVector& b);
};

#endif
