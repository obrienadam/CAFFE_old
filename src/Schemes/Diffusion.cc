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
                if(j < uJ)
                {
                    alpha = FvScheme::getAlpha(i, j, k, NORTH);

                    phiBar = alpha*gradPhiField_(i, j, k) + (1. - alpha)*gradPhiField_(i, j + 1, k);
                    A = mesh.fAreaNormN(i, j, k);
                    es = mesh.rnCellN(i, j, k);
                    phi1 = phiField(i, j + 1, k);
                    ds = mesh.rCellMagN(i, j, k);

                    gradPhiField_.faceN(i, j, k) = (phi1 - phi0)*A/(ds*dot(A, es)) + phiBar - dot(phiBar, es)*A/dot(A, es);
                }

                // Top faces
                if(k < uK)
                {
                    alpha = FvScheme::getAlpha(i, j, k, TOP);

                    phiBar = alpha*gradPhiField_(i, j, k) + (1. - alpha)*gradPhiField_(i, j, k + 1);
                    A = mesh.fAreaNormT(i, j, k);
                    es = mesh.rnCellT(i, j, k);
                    phi1 = phiField(i, j, k + 1);
                    ds = mesh.rCellMagT(i, j, k);

                    gradPhiField_.faceT(i, j, k) = (phi1 - phi0)*A/(ds*dot(A, es)) + phiBar - dot(phiBar, es)*A/dot(A, es);
                }
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
    for(k = 0; k < nCellsK; ++k)
    {
        for(i = 0; i < nCellsI; ++i)
        {
            phi0 = phiField(i, uJ, k);
            phiB = phiField(i, uJ + 1, k);
            db = mesh.rFaceMagN(i, uJ, k);
            eb = mesh.rnFaceN(i, uJ, k);
            A = mesh.fAreaNormN(i, uJ, k);

            gradPhiField_.faceN(i, uJ, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField_(i, uJ, k) - dot(gradPhiField_(i, uJ, k), eb)*A/dot(A, eb);

            phi0 = phiField(i, 0, k);
            phiB = phiField(i, -1, k);
            db = mesh.rFaceMagS(i, 0, k);
            eb = mesh.rnFaceS(i, 0, k);
            A = mesh.fAreaNormS(i, 0, k);

            gradPhiField_.faceS(i, 0, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField_(i, 0, k) - dot(gradPhiField_(i, 0, k), eb)*A/dot(A, eb);
        }
    }

    // Top and bottom
    for(j = 0; j < nCellsJ; ++j)
    {
        for(i = 0; i < nCellsI; ++i)
        {
            phi0 = phiField(i, j, uK);
            phiB = phiField(i, j, uK + 1);
            db = mesh.rFaceMagT(i, j, uK);
            eb = mesh.rnFaceT(i, j, uK);
            A = mesh.fAreaNormT(i, j, uK);

            gradPhiField_.faceT(i, j, uK) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField_(i, j, uK) - dot(gradPhiField_(i, j, uK), eb)*A/dot(A, eb);

            phi0 = phiField(i, j, 0);
            phiB = phiField(i, j, -1);
            db = mesh.rFaceMagB(i, j, 0);
            eb = mesh.rnFaceB(i, j, 0);
            A = mesh.fAreaNormB(i, j, 0);

            gradPhiField_.faceB(i, j, 0) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField_(i, j, 0) - dot(gradPhiField_(i, j, 0), eb)*A/dot(A, eb);
        }
    }
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

                if(i < nFacesI - 1 && k < nFacesK - 1)
                    phiField.fluxJ(i, j, k) = dot(gradPhiField_.faceJ(i, j, k), mesh.fAreaNormJ(i, j, k));

                if(i < nFacesI - 1 && j < nFacesJ - 1)
                    phiField.fluxK(i, j, k) = dot(gradPhiField_.faceK(i, j, k), mesh.fAreaNormK(i, j, k));
            }
        }
    }
}

// ************* Public Methods *************

void Diffusion::initialize(Input &input, HexaFvmMesh &mesh, std::string conservedFieldName)
{
    int nCellsI, nCellsJ, nCellsK;
    Matrix lsMatrix(6, 3);

    FvScheme::initialize(input, mesh, conservedFieldName);
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

void Diffusion::copySolution(std::vector<double> &original)
{
    int k, size;
    Field<double>& phiField = *phiFieldPtr_;

    size = phiField.size();

    for(k = 0; k < size; ++k)
    {
        original[k] = phiField(k);
    }
}

void Diffusion::updateSolution(std::vector<double>& update, int method)
{
    int k, size;
    Field<double>& phiField = *phiFieldPtr_;

    size = phiField.size();

    switch(method)
    {
    case ADD:

        for(k = 0; k < size; ++k)
        {
            phiField(k) += update[k];
        }

        break;
    case REPLACE:

        for(k = 0; k < size; ++k)
        {
            phiField(k) = update[k];
        }

        break;
    default:

        Output::raiseException("Diffusion", "updateSolution", "Invalid update method selected.");
    };
}
