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

// ************* Constructors and Destructors *************


Diffusion::Diffusion()
    :
      gradPhiField_("gradPhiField", PRIMITIVE)
{

}

Diffusion::~Diffusion()
{

}

// ************* Private Methods *************

void Diffusion::computeCellCenteredGradients()
{
    int i, j, k, nCellsI, nCellsJ, nCellsK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& phiField = *phiFieldPtr_;
    Matrix A(6, 3), b(6, 1);

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                b.reallocate(6, 1);

                A.addVector3DToRow(mesh.rCellE(i, j, k), 0, 0);
                A.addVector3DToRow(mesh.rCellW(i, j, k), 1, 0);
                A.addVector3DToRow(mesh.rCellN(i, j, k), 2, 0);
                A.addVector3DToRow(mesh.rCellS(i, j, k), 3, 0);
                A.addVector3DToRow(mesh.rCellT(i, j, k), 4, 0);
                A.addVector3DToRow(mesh.rCellB(i, j, k), 5, 0);

                b(0, 0) = phiField(i + 1, j, k) - phiField(i, j, k);
                b(1, 0) = phiField(i - 1, j, k) - phiField(i, j, k);
                b(2, 0) = phiField(i, j + 1, k) - phiField(i, j, k);
                b(3, 0) = phiField(i, j - 1, k) - phiField(i, j, k);
                b(4, 0) = phiField(i, j, k + 1) - phiField(i, j, k);
                b(5, 0) = phiField(i, j, k - 1) - phiField(i, j, k);

                A.solveLeastSquares(b);

                gradPhiField_(i, j, k).x = b(0, 0);
                gradPhiField_(i, j, k).y = b(1, 0);
                gradPhiField_(i, j, k).z = b(2, 0);
            }
        }
    }
}

void Diffusion::computeFaceCenteredGradients()
{
    int i, j, k, nCellsI, nCellsJ, nCellsK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& phiField = *phiFieldPtr_;
    Vector3D phiBar;
    double alpha;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    //- Reconstruct the interior faces

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                // East faces
                alpha = mesh.cellVol(i, j, k)/(mesh.cellVol(i, j, k) + mesh.cellVol(i + 1, j, k));
                phiBar = alpha*gradPhiField_(i, j, k) + (1. - alpha)*gradPhiField_(i + 1, j, k);
                gradPhiField_.faceE(i, j, k) = phiBar;

                // North faces
                alpha = mesh.cellVol(i, j, k)/(mesh.cellVol(i, j, k) + mesh.cellVol(i, j + 1, k));
                phiBar = alpha*gradPhiField_(i, j, k) + (1. - alpha)*gradPhiField_(i, j + 1, k);
                gradPhiField_.faceN(i, j, k) = phiBar;

                // Top faces
                alpha = mesh.cellVol(i, j, k)/(mesh.cellVol(i, j, k) + mesh.cellVol(i, j, k + 1));
                phiBar = alpha*gradPhiField_(i, j, k) + (1. - alpha)*gradPhiField_(i, j, k + 1);
                gradPhiField_.facesT(i, j, k) = phiBar;
            }
        }
    }
}

// ************* Public Methods *************

void Diffusion::initialize(HexaFvmMesh &mesh, std::string conservedFieldName)
{
    int nCellsI, nCellsJ, nCellsK;
    Matrix lsMatrix(6, 3);

    FvScheme::initialize(mesh, conservedFieldName);
    phiFieldPtr_ = &mesh.findScalarField(conservedFieldName_);

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    gradPhiField_.allocate(nCellsI, nCellsJ, nCellsK);
    stencil_.allocate(3, 3, 3);
}

int Diffusion::nConservedVariables()
{
    return phiFieldPtr_->size();
}

void Diffusion::discretize(std::vector<double>& timeDerivatives)
{
    int i, j, k, l, nCellsI, nCellsJ, nCellsK;
    Field<double>& phiField = *phiFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;
    Vector3D gradPhiBar, gradPhiFace;

    nCellsI = mesh.nCellsI() - 1;
    nCellsJ = mesh.nCellsJ() - 1;
    nCellsK = mesh.nCellsK() - 1;

    computeCellCenteredGradients();
    computeFaceCenteredGradients();
    //gradPhiField_.print();

    for(k = 0, l = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i, ++l)
            {

                timeDerivatives[l] = phiField.sumFluxes(i, j, k)/mesh.cellVol(i, j, k);
            }
        }
    }

}

void Diffusion::updateSolution(std::vector<double>& update)
{
    int i, j, k, l, nCellsI, nCellsJ, nCellsK;
    Field<double>& phiField = *phiFieldPtr_;

    nCellsI = phiField.sizeI();
    nCellsJ = phiField.sizeJ();
    nCellsK = phiField.sizeK();

    for(k = 0, l = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i, ++l)
            {
                phiField(i, j, k) = update[l];
            }
        }
    }
}
