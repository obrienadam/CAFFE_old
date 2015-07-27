/**
 * @file    Diffusion.cc
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
 * This file contains the implementations for the methods of class
 * Diffusion.
 */

#include <math.h>

#include "Diffusion.h"
#include "Output.h"
#include "FvScalarScheme.h"

// ************* Constructors and Destructors *************


Diffusion::Diffusion(const Input &input, const HexaFvmMesh &mesh)
    :
      Solver(input, mesh),
      phiField_(mesh, Field<double>::CONSERVED, "phi"),
      gradPhiField_(mesh, Field<Vector3D>::AUXILLARY, "gradPhi"),
      bcs_(input, phiField_)
{
    using namespace std;

    Solver::createMatrices(1, 1, 7);
    mesh_.addArray3DToTecplotOutput(phiField_.name, phiField_.cellData());
    mesh_.addArray3DToTecplotOutput(gradPhiField_.name, gradPhiField_.cellData());
}

Diffusion::~Diffusion()
{

}

// ************* Public Methods *************

double Diffusion::solve(double time_Step)
{
    int i, j, k, cols[7];
    double a0P, a[7], b;

    time_.tic();
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(Solver::solutionType_ == Solver::UNSTEADY)
                    a0P = 1./time_Step;
                else
                    a0P = 0.;

                a[1] = mesh_.dE(i, j, k);
                a[2] = mesh_.dW(i, j, k);
                a[3] = mesh_.dN(i, j, k);
                a[4] = mesh_.dS(i, j, k);
                a[5] = mesh_.dT(i, j, k);
                a[6] = mesh_.dB(i, j, k);

                a[0] = -(a[1] + a[2] + a[3] + a[4] + a[5] + a[6]) + a0P;
                b = a0P*phiField_(i, j, k);

                bcs_.setImplicitBoundaryCoefficients(i, j, k, a, b);

                cols[0] = indexMap_(i, j, k, 0);
                cols[1] = indexMap_(i + 1, j, k, 0);
                cols[2] = indexMap_(i - 1, j, k, 0);
                cols[3] = indexMap_(i, j + 1, k, 0);
                cols[4] = indexMap_(i, j - 1, k, 0);
                cols[5] = indexMap_(i, j, k + 1, 0);
                cols[6] = indexMap_(i, j, k - 1, 0);

                A_[0].setRow(indexMap_(i, j, k, 0), 7, cols, a);
                b_[0].setValue(indexMap_(i, j, k, 0), b);
            }
        }
    }

    time_.toc();
    Output::print("\nDiffusion", "matrix assembly completed in " + time_.elapsedTime());

    biCGStabIters_ = A_[0].solve(b_[0], x_[0]);

    time_.tic();
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                phiField_(i, j, k) = x_[0](indexMap_(i, j, k, 0));
            }
        }
    }

    time_.toc();
    Output::print("Diffusion", "solution of linear system completed in " + time_.elapsedTime());

    bcs_.setBoundaries();
    FvScalarScheme::extrapolateInteriorFaces(FvScheme::DIVERGENCE_THEOREM, phiField_, gradPhiField_);
    FvScalarScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, phiField_, gradPhiField_);

    scale(-1., b_[0]);
    multiplyAdd(A_[0], x_[0], b_[0], res_[0]);

    return res_[0].infNorm();
}
