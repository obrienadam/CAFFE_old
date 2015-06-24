/**
 * @file    SparseMatrix.cc
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

SparseMatrix::SparseMatrix()
    :
      maxIters_(50000),
      rToler_(1e-5),
      absToler_(1e-8)
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
    MatDestroy(&A_);
    KSPDestroy(&ksp_);
}

void SparseMatrix::allocate(int m, int n, int nnz)
{
    MatCreate(PETSC_COMM_WORLD, &A_);
    MatSetSizes(A_, PETSC_DECIDE, PETSC_DECIDE, m, m);
    // Should be MATMPIAIJ if parallel or MATSEQAIJ if single
    MatSetType(A_, MATAIJ);
    // Both of these need to be called to ensure it works with MPI and single processor (very weird)
    MatMPIAIJSetPreallocation(A_, nnz, NULL, nnz, NULL);
    MatSeqAIJSetPreallocation(A_, nnz, NULL);
    MatSetUp(A_);
}

void SparseMatrix::setValue(int i, int j, double value, InsertMode insertMode)
{
    MatSetValues(A_, 1, &i, 1, &j, &value, insertMode);
}

int SparseMatrix::solve(const SparseVector& b, SparseVector& x)
{
    int nIters;

    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

    KSPCreate(PETSC_COMM_WORLD, &ksp_);
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

void SparseMatrix::print()
{
    MatView(A_, PETSC_VIEWER_STDOUT_SELF);
}
