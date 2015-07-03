/**
 * @file    SparseMatrix.h
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
 * Contains the interface for class SparseMatrix, which is for sparse
 * storage of large linear systems. This class is built using the Petsc
 * library.
 */

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <petsc.h>

#include "SparseVector.h"
#include "Input.h"

class SparseMatrix
{
private:

    Mat A_;
    KSP ksp_;
    PC pc_;

    int maxIters_;
    double rToler_, absToler_;
    PetscErrorCode errorCode_;

    // MPI related
    int iLower_, iUpper_;

public:

    SparseMatrix();
    SparseMatrix(int m, int n, int nnz);
    ~SparseMatrix();

    void allocate(int m, int n, int nnz);
    void deallocate();
    void setValue(int i, int j, double value);
    void setRow(int rowNo, int nCols, int colNos[], double values[]);
    void zeroEntries();

    int solve(const SparseVector &b, SparseVector& x);
    void assemble();

    void print();

    friend void multiply(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &C);
    friend void multiply(const SparseMatrix &A, const SparseVector &x, SparseVector &b);
    friend void multiplyAdd(const SparseMatrix& A, const SparseVector& x1, const SparseVector& x2, SparseVector& b);
    friend void scale(double alpha, SparseMatrix& A);
};

void multiply(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& C);
void multiply(const SparseMatrix& A, const SparseVector& x, SparseVector& b);
void multiplyAdd(const SparseMatrix& A, const SparseVector& x1, const SparseVector& x2, SparseVector& b);
void scale(double alpha, SparseMatrix& A);

#endif
