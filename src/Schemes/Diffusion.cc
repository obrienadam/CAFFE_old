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

#include "Diffusion.h"
#include "Output.h"

// ************* Private Methods *************

void Diffusion::computeCellCenteredGradients()
{
    int i, j, k, nCellsI, nCellsJ, nCellsK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& phiField = *phiFieldPtr_;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                bLs_(0, 0) = phiField(i + 1, j, k) - phiField(i, j, k);
                bLs_(1, 0) = phiField(i - 1, j, k) - phiField(i, j, k);
                bLs_(2, 0) = phiField(i, j + 1, k) - phiField(i, j, k);
                bLs_(3, 0) = phiField(i, j - 1, k) - phiField(i, j, k);
                bLs_(4, 0) = phiField(i, j, k + 1) - phiField(i, j, k);
                bLs_(5, 0) = phiField(i, j, k - 1) - phiField(i, j, k);

                xLs_ = solveLeastSquares(lsMatrices_(i, j, k), bLs_);

                gradPhiField_(i, j, k).x = xLs_(0, 0);
                gradPhiField_(i, j, k).y = xLs_(1, 0);
                gradPhiField_(i, j, k).z = xLs_(2, 0);
            }
        }
    }
}

// ************* Public Methods *************

Diffusion::Diffusion()
{

}

Diffusion::~Diffusion()
{

}

void Diffusion::initialize(HexaFvmMesh &mesh, std::string conservedFieldName)
{
    int i, j, k, nCellsI, nCellsJ, nCellsK;
    Vector3D tmp;
    Matrix lsMatrix(6, 3);

    FvScheme::initialize(mesh, conservedFieldName);
    phiFieldPtr_ = &mesh.findScalarField(conservedFieldName_);

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    lsMatrices_.allocate(nCellsI, nCellsJ, nCellsK);

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                lsMatrix.addVector3DToRow(mesh.nesE(i, j, k)*mesh.cellToCellDistanceE(i, j, k), 0, 0);
                lsMatrix.addVector3DToRow(mesh.nesW(i, j, k)*mesh.cellToCellDistanceW(i, j, k), 1, 0);
                lsMatrix.addVector3DToRow(mesh.nesN(i, j, k)*mesh.cellToCellDistanceN(i, j, k), 2, 0);
                lsMatrix.addVector3DToRow(mesh.nesS(i, j, k)*mesh.cellToCellDistanceS(i, j, k), 3, 0);
                lsMatrix.addVector3DToRow(mesh.nesT(i, j, k)*mesh.cellToCellDistanceT(i, j, k), 4, 0);
                lsMatrix.addVector3DToRow(mesh.nesB(i, j, k)*mesh.cellToCellDistanceB(i, j, k), 5, 0);

                lsMatrices_(i, j, k) = lsMatrix;
            }
        }
    }

    aLs_.allocate(6, 3);
    bLs_.allocate(6, 1);
    aLs_.allocate(3, 1);
    gradPhiField_.allocate(nCellsI, nCellsJ, nCellsK);
    stencil_.allocate(3, 3, 3);
}

int Diffusion::nConservedVariables()
{
    return phiFieldPtr_->size();
}

void Diffusion::discretize(std::vector<double>& timeDerivatives)
{
    Field<double>& phiField = *phiFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;

    computeCellCenteredGradients();

}

void Diffusion::updateSolution(std::vector<double> &timeDerivatives)
{

}
