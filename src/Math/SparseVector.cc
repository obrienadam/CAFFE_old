#include <iostream>

#include "SparseVector.h"

SparseVector::SparseVector()
    :
      m_(0)
{
    PetscInitializeNoArguments();
    errorCode_ = VecCreate(PETSC_COMM_WORLD, &vec_);
    errorCode_ = VecSetType(vec_, VECMPI);
}

SparseVector::SparseVector(int m)
    :
      SparseVector()
{
    setSize(m);
}

SparseVector::~SparseVector()
{
    VecDestroy(&vec_);
}

void SparseVector::setSize(int m)
{
    m_ = m;
    VecSetSizes(vec_, PETSC_DECIDE, m_);
    VecGetOwnershipRange(vec_, &iStart_, &iEnd_);
}

void SparseVector::setValue(int i, double value, InsertMode insertMode)
{
    VecSetValue(vec_, i, value, insertMode);
}

void SparseVector::setValues(int n, int *indices, double *values, InsertMode insertMode)
{
    VecSetValues(vec_, n, indices, values, insertMode);
}

double SparseVector::operator ()(int i)
{
    double value;

    errorCode_ = VecGetValues(vec_, 1, &i, &value);

    return value;
}

void SparseVector::beginAssembly()
{
    VecAssemblyBegin(vec_);
}

void SparseVector::endAssembly()
{
    VecAssemblyEnd(vec_);
}

void SparseVector::assemble()
{
    VecAssemblyBegin(vec_);
    VecAssemblyEnd(vec_);
}

void SparseVector::print()
{
    PetscScalar *array;

    VecGetArray(vec_, &array);

    for(int i = iStart_; i < iEnd_; ++i)
    {
        std::cout << array[i] << std::endl;
    }

    VecRestoreArray(vec_, &array);
}
