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
      rToler_(1e-8),
      absToler_(1e-8)
{
    PetscInitializeNoArguments();
    errorCode_ = MatCreate(PETSC_COMM_WORLD, &A_);
    errorCode_ = MatSetType(A_, MATMPIAIJ);
    KSPCreate(PETSC_COMM_WORLD, &ksp_);
}

SparseMatrix::SparseMatrix(int m, int n)
    :
      SparseMatrix()
{
    setSize(m, n);
}

SparseMatrix::~SparseMatrix()
{
    MatDestroy(&A_);
    KSPDestroy(&ksp_);
}

void SparseMatrix::initialize(Input &input)
{
    maxIters_ = input.inputInts["kspMaxItrs"];
    rToler_ = input.inputDoubles["kspRelativeTolerance"];
    absToler_ = input.inputDoubles["kspAbsoluteTolerance"];
}

void SparseMatrix::setSize(int m, int n)
{
    m_ = m;
    n_ = n;
    MatSetSizes(A_, PETSC_DECIDE, PETSC_DECIDE, m_, n_);
    errorCode_ = MatSetUp(A_);
}

void SparseMatrix::preallocate(int dnz, int onz)
{
    errorCode_ = MatMPIAIJSetPreallocation(A_, dnz, NULL, onz, NULL);
}

void SparseMatrix::setValue(int i, int j, double value, InsertMode insertMode)
{
    MatSetValue(A_, i, j, value, insertMode);
}

void SparseMatrix::setValues(int m, int *iIndices, int n, int *jIndices, double *values, InsertMode insertMode)
{
    MatSetValues(A_, m, iIndices, n, jIndices, values, insertMode);
}

void SparseMatrix::beginAssembly()
{
    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
}

void SparseMatrix::endAssembly()
{
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
    KSPSetOperators(ksp_, A_, A_);
    KSPGetPC(ksp_, &pc_);
    PCSetType(pc_, PCJACOBI);
    KSPSetTolerances(ksp_, rToler_, absToler_, PETSC_DEFAULT, maxIters_);
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
    KSPSetUp(ksp_);
}

void SparseMatrix::assemble()
{
    beginAssembly();
    endAssembly();
}

int SparseMatrix::solve(const SparseVector& b, SparseVector& x)
{
    int nIters;

    KSPSolve(ksp_, b.vec_, x.vec_);
    KSPGetIterationNumber(ksp_, &nIters);

    return nIters;
}

void SparseMatrix::print()
{
    MatView(A_, PETSC_VIEWER_STDOUT_SELF);
}
