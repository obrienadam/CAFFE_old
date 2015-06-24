#include <iostream>

#include "SparseVector.h"

SparseVector::SparseVector()
{
    PetscInitializeNoArguments();
}

SparseVector::SparseVector(int m)
    :
      SparseVector()
{
    allocate(m);
}

SparseVector::~SparseVector()
{
    VecDestroy(&vec_);
}

void SparseVector::allocate(int m)
{
    VecCreate(PETSC_COMM_WORLD, &vec_);
    VecSetSizes(vec_, PETSC_DECIDE, m);
    VecSetType(vec_, VECSTANDARD);
}

void SparseVector::allocate(const SparseVector& other)
{
    VecDuplicate(other.vec_, &vec_);
}

void SparseVector::setValue(int i, double value, InsertMode insertMode)
{
    VecSetValues(vec_, 1, &i, &value, insertMode);
}

double SparseVector::operator ()(int i)
{
    double value;

    VecGetValues(vec_, 1, &i, &value);

    return value;
}

void SparseVector::print()
{
    VecView(vec_, PETSC_VIEWER_STDOUT_SELF);
}
