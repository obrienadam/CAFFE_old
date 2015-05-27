#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H

#include <petsc.h>

class SparseVector
{
private:

    int m_;
    int iStart_, iEnd_;

    Vec vec_;

    PetscErrorCode errorCode_;

public:

    SparseVector();
    SparseVector(int m);
    ~SparseVector();

    void setSize(int m);

    void setValue(int i, double value, InsertMode insertMode = INSERT_VALUES);
    void setValues(int n, int* indices, double* values, InsertMode insertMode = INSERT_VALUES);

    double operator()(int i);

    void beginAssembly();
    void endAssembly();
    void assemble();

    void print();

    friend class SparseMatrix;
};

#endif
