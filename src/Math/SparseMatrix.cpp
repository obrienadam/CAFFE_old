/**
 * @file    SparseMatrix.cpp
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
 * Contains the implementations of methods for class SparseMatrix.
 */

#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(double rToler, double absToler, int maxIters)
    :
      rToler_(rToler),
      absToler_(absToler),
      maxIters_(maxIters)
{
    PetscInitializeNoArguments();
}

SparseMatrix::SparseMatrix(int m, int n, int nnz)
    :
      SparseMatrix()
{
    allocate(m, n, nnz);
}

SparseMatrix::~SparseMatrix()
{
    deallocate();
}

void SparseMatrix::allocate(int m, int n, int nnz)
{
    MatCreate(PETSC_COMM_WORLD, &A_);
    MatSetSizes(A_, PETSC_DECIDE, PETSC_DECIDE, m, n);
    // Should be MATMPIAIJ if parallel or MATSEQAIJ if single
    MatSetType(A_, MATAIJ);
    // Both of these need to be called to ensure it works with MPI and single processor (very weird)
    MatMPIAIJSetPreallocation(A_, nnz, NULL, 0, NULL);
    MatSeqAIJSetPreallocation(A_, nnz, NULL);
    MatSetUp(A_);

    // Get MPI ranges
    MatGetOwnershipRange(A_, &iLower_, &iUpper_);

    // Create solver context
    KSPCreate(PETSC_COMM_WORLD, &ksp_);
}

void SparseMatrix::deallocate()
{
    MatDestroy(&A_);
    KSPDestroy(&ksp_);
}

void SparseMatrix::setValue(int i, int j, double value)
{
    MatSetValues(A_, 1, &i, 1, &j, &value, INSERT_VALUES);
}

void SparseMatrix::setRow(int rowNo, int nCols, int colNos[], double values[])
{
    MatSetValues(A_, 1, &rowNo, nCols, colNos, values, INSERT_VALUES);
}

void SparseMatrix::zeroEntries()
{
    MatZeroEntries(A_);
}

int SparseMatrix::solve(const SparseVector& b, SparseVector& x)
{
    int nIters;

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
    KSPSetOperators(ksp_, A_, A_);
    KSPGetPC(ksp_, &pc_);
    PCSetType(pc_, PCILU);
    PCFactorSetFill(pc_, 2);
    KSPSetTolerances(ksp_, rToler_, absToler_, PETSC_DEFAULT, maxIters_);
    KSPSetType(ksp_, KSPBCGS);

    KSPSolve(ksp_, b.vec_, x.vec_);
    KSPGetIterationNumber(ksp_, &nIters);

    return nIters;
}

void SparseMatrix::assemble()
{
    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
}

void SparseMatrix::print()
{
    MatView(A_, PETSC_VIEWER_STDOUT_SELF);
}

void multiply(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &C)
{
    MatMatMult(A.A_, B.A_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C.A_);
}

void multiply(const SparseMatrix &A, const SparseVector &x, SparseVector &b)
{
    MatMult(A.A_, x.vec_, b.vec_);
}

void multiplyAdd(const SparseMatrix &A, const SparseVector &x1, const SparseVector &x2, SparseVector &b)
{
    MatMultAdd(A.A_, x1.vec_, x2.vec_, b.vec_);
}

void scale(double alpha, SparseMatrix &A)
{
    MatScale(A.A_, alpha);
}
