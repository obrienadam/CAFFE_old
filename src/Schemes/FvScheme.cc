/**
 * @file    FvScheme.cc
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
 * This file contains the implementations of methods for class FvScheme.
 */

#include <cstdlib>
#include <math.h>

#include "FvScheme.h"
#include "Output.h"
#include "Matrix.h"

FvScheme::FvScheme()
    :
      conservedFieldName_("phi"),
      meshPtr_(NULL)
{

}

void FvScheme::initialize(Input &input, HexaFvmMesh &mesh, std::string conservedFieldName)
{
    meshPtr_ = &mesh;
    nCellsI_ = mesh.nCellsI();
    nCellsJ_ = mesh.nCellsJ();
    nCellsK_ = mesh.nCellsK();
    nFacesI_ = mesh.nFacesI();
    nFacesJ_ = mesh.nFacesJ();
    nFacesK_ = mesh.nFacesK();

    conservedFieldName_ = conservedFieldName;
}

double FvScheme::getAlpha(int i, int j, int k, int direction)
{
    int deltaI = 0, deltaJ = 0, deltaK = 0;

    switch (direction)
    {
    case EAST:
        ++deltaI;
        break;
    case WEST:
        --deltaI;
        break;
    case NORTH:
        ++deltaJ;
        break;
    case SOUTH:
        --deltaJ;
        break;
    case TOP:
        ++deltaK;
        break;
    case BOTTOM:
        --deltaK;
        break;
    default:
        Output::raiseException("FvScheme", "getAlpha", "Invalid direction specified.");
    };

    //- This is a temporary fix and should probably be looked at again in the future.

    if(i + deltaI > 0 && i + deltaI < nCellsI_ - 1
            && j + deltaJ > 0 && j + deltaJ < nCellsJ_ - 1
            && k + deltaK > 0 && k + deltaK < nCellsK_ - 1)
        return meshPtr_->cellVol(i, j, k)/(meshPtr_->cellVol(i, j, k) + meshPtr_->cellVol(i + deltaI, j + deltaJ, k + deltaK));
    else
        return 1.;
}

void FvScheme::computeCellCenteredGradients(Field<double> &phiField, Field<Vector3D> &gradPhiField, int method)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;
    Matrix A(6, 3), b(6, 1);

    switch(method)
    {
    case LEAST_SQUARES:

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
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

                gradPhiField(i, j, k).x = b(0, 0);
                gradPhiField(i, j, k).y = b(1, 0);
                gradPhiField(i, j, k).z = b(2, 0);
            }
        }
    }

        break;
    case DIVERGENCE_THEOREM:

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    gradPhiField(i, j, k) = 1./mesh.cellVol(i, j, k)*(0.5*(phiField(i, j, k) + phiField(i + 1, j, k))*mesh.fAreaNormE(i, j, k)
                                                                      + 0.5*(phiField(i, j, k) + phiField(i - 1, j, k))*mesh.fAreaNormW(i, j, k)
                                                                      + 0.5*(phiField(i, j, k) + phiField(i, j + 1, k))*mesh.fAreaNormN(i, j, k)
                                                                      + 0.5*(phiField(i, j, k) + phiField(i, j - 1, k))*mesh.fAreaNormS(i, j, k)
                                                                      + 0.5*(phiField(i, j, k) + phiField(i, j, k + 1))*mesh.fAreaNormT(i, j, k)
                                                                      + 0.5*(phiField(i, j, k) + phiField(i, j, k - 1))*mesh.fAreaNormB(i, j, k));
                }
            }
        }

        break;
    };
}

void FvScheme::computeCellCenteredJacobians(Field<Vector3D> &vecField, Field<Tensor3D> &tensorField, int method)
{
    int i, j, k, l, m;
    HexaFvmMesh& mesh = *meshPtr_;
    Matrix A(6, 3), b(6, 1), x(3, 1);

    switch(method)
    {
    case LEAST_SQUARES:

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                b.reallocate(6, 1);

                //- Construct the least-squares coefficient matrix for the cell

                A.addVector3DToRow(mesh.rCellE(i, j, k), 0, 0);
                A.addVector3DToRow(mesh.rCellW(i, j, k), 1, 0);
                A.addVector3DToRow(mesh.rCellN(i, j, k), 2, 0);
                A.addVector3DToRow(mesh.rCellS(i, j, k), 3, 0);
                A.addVector3DToRow(mesh.rCellT(i, j, k), 4, 0);
                A.addVector3DToRow(mesh.rCellB(i, j, k), 5, 0);

                for(l = 0; l < 3; ++l)
                {
                    b(0, 0) = vecField(i + 1, j, k)(l) - vecField(i, j, k)(l);
                    b(1, 0) = vecField(i - 1, j, k)(l) - vecField(i, j, k)(l);
                    b(2, 0) = vecField(i, j + 1, k)(l) - vecField(i, j, k)(l);
                    b(3, 0) = vecField(i, j - 1, k)(l) - vecField(i, j, k)(l);
                    b(4, 0) = vecField(i, j, k + 1)(l) - vecField(i, j, k)(l);
                    b(5, 0) = vecField(i, j, k - 1)(l) - vecField(i, j, k)(l);

                    x = solveLeastSquares(A, b);

                    for(m = 0; m < 3; ++m)
                    {
                        tensorField(i, j, k)(l, m) = x(m, 0);
                    }
                }

            }
        }
    }

        break;
    case DIVERGENCE_THEOREM:

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    tensorField(i, j, k) = 1./mesh.cellVol(i, j, k)*(tensor(0.5*(vecField(i, j, k) + vecField(i + 1, j, k)), mesh.fAreaNormE(i, j, k))
                                                                     + tensor(0.5*(vecField(i, j, k) + vecField(i - 1, j, k)), mesh.fAreaNormW(i, j, k))
                                                                     + tensor(0.5*(vecField(i, j, k) + vecField(i, j + 1, k)), mesh.fAreaNormN(i, j, k))
                                                                     + tensor(0.5*(vecField(i, j, k) + vecField(i, j - 1, k)), mesh.fAreaNormS(i, j, k))
                                                                     + tensor(0.5*(vecField(i, j, k) + vecField(i, j, k + 1)), mesh.fAreaNormT(i, j, k))
                                                                     + tensor(0.5*(vecField(i, j, k) + vecField(i, j, k - 1)), mesh.fAreaNormB(i, j, k)));
                }
            }
        }

        break;
    };
}

void FvScheme::computeFaceCenteredGradients(Field<double> &phiField, Field<Vector3D> &gradPhiField)
{
    int i, j, k, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;
    Vector3D phiBar, A, es, eb;
    double alpha, phi1, phi0, phiB, ds, db;

    uI = nCellsI_ - 1;
    uJ = nCellsJ_ - 1;
    uK = nCellsK_ - 1;

    //- Reconstruct the interior faces

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                phi0 = phiField(i, j, k);

                // East faces
                if(i < uI)
                {
                    alpha = FvScheme::getAlpha(i, j, k, EAST);

                    phiBar = alpha*gradPhiField(i, j, k) + (1. - alpha)*gradPhiField(i + 1, j, k);
                    A = mesh.fAreaNormE(i, j, k);
                    es = mesh.rnCellE(i, j, k);
                    phi1 = phiField(i + 1, j, k);
                    ds = mesh.rCellMagE(i, j, k);

                    gradPhiField.faceE(i, j, k) = (phi1 - phi0)*A/(ds*dot(A, es)) + phiBar - dot(phiBar, es)*A/dot(A, es);
                }

                // North faces
                if(j < uJ)
                {
                    alpha = FvScheme::getAlpha(i, j, k, NORTH);

                    phiBar = alpha*gradPhiField(i, j, k) + (1. - alpha)*gradPhiField(i, j + 1, k);
                    A = mesh.fAreaNormN(i, j, k);
                    es = mesh.rnCellN(i, j, k);
                    phi1 = phiField(i, j + 1, k);
                    ds = mesh.rCellMagN(i, j, k);

                    gradPhiField.faceN(i, j, k) = (phi1 - phi0)*A/(ds*dot(A, es)) + phiBar - dot(phiBar, es)*A/dot(A, es);
                }

                // Top faces
                if(k < uK)
                {
                    alpha = FvScheme::getAlpha(i, j, k, TOP);

                    phiBar = alpha*gradPhiField(i, j, k) + (1. - alpha)*gradPhiField(i, j, k + 1);
                    A = mesh.fAreaNormT(i, j, k);
                    es = mesh.rnCellT(i, j, k);
                    phi1 = phiField(i, j, k + 1);
                    ds = mesh.rCellMagT(i, j, k);

                    gradPhiField.faceT(i, j, k) = (phi1 - phi0)*A/(ds*dot(A, es)) + phiBar - dot(phiBar, es)*A/dot(A, es);
                }
            }
        }
    }

    //- Reconstruct the boundary faces

    // East and west
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            phi0 = phiField(uI, j, k);
            phiB = phiField(uI + 1, j, k);
            db = mesh.rFaceMagE(uI, j, k);
            eb = mesh.rnFaceE(uI, j, k);
            A = mesh.fAreaNormE(uI, j, k);

            gradPhiField.faceE(uI, j, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField(uI, j, k) - dot(gradPhiField(uI, j, k), eb)*A/dot(A, eb);

            phi0 = phiField(0, j, k);
            phiB = phiField(-1, j, k);
            db = mesh.rFaceMagW(0, j, k);
            eb = mesh.rnFaceW(0, j, k);
            A = mesh.fAreaNormW(0, j, k);

            gradPhiField.faceW(0, j, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField(0, j, k) - dot(gradPhiField(0, j, k), eb)*A/dot(A, eb);
        }
    }

    // North and south
    for(k = 0; k < nCellsK_; ++k)
    {
        for(i = 0; i < nCellsI_; ++i)
        {
            phi0 = phiField(i, uJ, k);
            phiB = phiField(i, uJ + 1, k);
            db = mesh.rFaceMagN(i, uJ, k);
            eb = mesh.rnFaceN(i, uJ, k);
            A = mesh.fAreaNormN(i, uJ, k);

            gradPhiField.faceN(i, uJ, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField(i, uJ, k) - dot(gradPhiField(i, uJ, k), eb)*A/dot(A, eb);

            phi0 = phiField(i, 0, k);
            phiB = phiField(i, -1, k);
            db = mesh.rFaceMagS(i, 0, k);
            eb = mesh.rnFaceS(i, 0, k);
            A = mesh.fAreaNormS(i, 0, k);

            gradPhiField.faceS(i, 0, k) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField(i, 0, k) - dot(gradPhiField(i, 0, k), eb)*A/dot(A, eb);
        }
    }

    // Top and bottom
    for(j = 0; j < nCellsJ_; ++j)
    {
        for(i = 0; i < nCellsI_; ++i)
        {
            phi0 = phiField(i, j, uK);
            phiB = phiField(i, j, uK + 1);
            db = mesh.rFaceMagT(i, j, uK);
            eb = mesh.rnFaceT(i, j, uK);
            A = mesh.fAreaNormT(i, j, uK);

            gradPhiField.faceT(i, j, uK) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField(i, j, uK) - dot(gradPhiField(i, j, uK), eb)*A/dot(A, eb);

            phi0 = phiField(i, j, 0);
            phiB = phiField(i, j, -1);
            db = mesh.rFaceMagB(i, j, 0);
            eb = mesh.rnFaceB(i, j, 0);
            A = mesh.fAreaNormB(i, j, 0);

            gradPhiField.faceB(i, j, 0) = (phiB - phi0)*A/(db*dot(A, eb)) + gradPhiField(i, j, 0) - dot(gradPhiField(i, j, 0), eb)*A/dot(A, eb);
        }
    }
}

void FvScheme::getMeshStencil(int i, int j, int k, int direction, Vector3D &faceNorm, Vector3D &cellRelVec, double& alpha)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int deltaI = 0, deltaJ = 0, deltaK = 0;

    switch (direction)
    {
    case EAST:

        faceNorm = mesh.fAreaNormE(i, j, k);
        cellRelVec = mesh.rCellE(i, j, k);
        ++deltaI;

        break;
    case WEST:

        faceNorm = mesh.fAreaNormW(i, j, k);
        cellRelVec = mesh.rCellW(i, j, k);
        --deltaI;

        break;
    case NORTH:

        faceNorm = mesh.fAreaNormN(i, j, k);
        cellRelVec = mesh.rCellN(i, j, k);
        ++deltaJ;

        break;
    case SOUTH:

        faceNorm = mesh.fAreaNormS(i, j, k);
        cellRelVec = mesh.rCellS(i, j, k);
        --deltaJ;

        break;
    case TOP:

        faceNorm = mesh.fAreaNormT(i, j, k);
        cellRelVec = mesh.rCellT(i, j, k);
        ++deltaK;

        break;

    case BOTTOM:

        faceNorm = mesh.fAreaNormB(i, j, k);
        cellRelVec = mesh.rCellB(i, j, k);
        --deltaK;

        break;

    default:

        Output::raiseException("FvScheme", "getMeshStencil", "Invalid direction specified.");
    };

    if(i + deltaI > 0 && i + deltaI < nCellsI_ - 1
            && j + deltaJ > 0 && j + deltaJ < nCellsJ_ - 1
            && k + deltaK > 0 && k + deltaK < nCellsK_ - 1)
    {
        alpha = mesh.cellVol(i, j, k)/(mesh.cellVol(i, j, k) + mesh.cellVol(i + deltaI, j + deltaJ, k + deltaK));
    }
    else
        alpha = 1.;
}

void FvScheme::displayUpdateMessage()
{

}
