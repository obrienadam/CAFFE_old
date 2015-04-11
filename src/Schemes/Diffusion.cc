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
    int i, j, k, nCellsI, nCellsJ, nCellsK, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& phiField = *phiFieldPtr_;
    Vector3D phiBar, A, es, eb;
    double alpha, phi1, phi0, phiB, ds, db;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    uI = nCellsI - 1;
    uJ = nCellsJ - 1;
    uK = nCellsK - 1;

    //- Reconstruct the interior faces

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                phi0 = phiField(i, j, k);

                // East faces
                if(i < uI)
                {
                    alpha = FvScheme::getAlpha(i, j, k, EAST);

                    phiBar = alpha*gradPhiField_(i, j, k) + (1. - alpha)*gradPhiField_(i + 1, j, k);
                    A = mesh.fAreaNormE(i, j, k);
                    es = mesh.rnCellE(i, j, k);
                    phi1 = phiField(i + 1, j, k);
                    ds = mesh.rCellMagE(i, j, k);

                    gradPhiField_.faceE(i, j, k) = (phi1 - phi0)*A/(ds*dot(A, es)) + phiBar - dot(phiBar, es)*A/dot(A, es);
                }

                // North faces

                // Top faces
            }
        }
    }

    //- Reconstruct the boundary faces

    // East and west
    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            phi0 = phiField(uI, j, k);
            phiB = phiField(uI + 1, j, k);
            db = mesh.rFaceMagE(uI, j, k);
            eb = mesh.rnFaceE(uI, j, k);
            A = mesh.fAreaNormE(uI, j, k);

            gradPhiField_.faceE(uI, j, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField_(uI, j, k) - dot(gradPhiField_(uI, j, k), eb)*A/dot(A, eb);

            phi0 = phiField(0, j, k);
            phiB = phiField(-1, j, k);
            db = mesh.rFaceMagW(0, j, k);
            eb = mesh.rnFaceW(0, j, k);
            A = mesh.fAreaNormW(0, j, k);

            gradPhiField_.faceW(0, j, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField_(0, j, k) - dot(gradPhiField_(0, j, k), eb)*A/dot(A, eb);
        }
    }

    // North and south

    // Top and bottom
}

void Diffusion::computeFaceFluxes()
{
    int i, j, k, nFacesI, nFacesJ, nFacesK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& phiField = *phiFieldPtr_;

    nFacesI = mesh.nFacesI();
    nFacesJ = mesh.nFacesJ();
    nFacesK = mesh.nFacesK();

    for(k = 0; k < nFacesK; ++k)
    {
        for(j = 0; j < nFacesJ; ++j)
        {
            for(i = 0; i < nFacesI; ++i)
            {
                if(j < nFacesJ - 1 && k < nFacesK - 1)
                    phiField.fluxI(i, j, k) = dot(gradPhiField_.faceI(i, j, k), mesh.fAreaNormI(i, j, k));

                //phiField.fluxJ(i, j, k) = dot(gradPhiField_.faceJ(i, j, k), mesh.fAreaNormJ(i, j, k));
                //phiField.fluxK(i, j, k) = dot(gradPhiField_.faceK(i, j, k), mesh.fAreaNormK(i, j, k));
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
    mesh.addVectorField("phiGrad");

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

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    phiField.setBoundaryFields();
    computeCellCenteredGradients();
    computeFaceCenteredGradients();
    computeFaceFluxes();

    for(k = 0, l = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                timeDerivatives[l] = phiField.sumFluxes(i, j, k)/mesh.cellVol(i, j, k);
                mesh.findVectorField("phiGrad")(i, j, k) = gradPhiField_(i, j, k);
                ++l;
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
            for(i = 0; i < nCellsI; ++i)
            {
                phiField(i, j, k) += update[l];
                ++l;
            }
        }
    }

    Output::print("Diffusion", "Update norm -> " + std::to_string(FvScheme::computeUpdateNorm(update)));
}
