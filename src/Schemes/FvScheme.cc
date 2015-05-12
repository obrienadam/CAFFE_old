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

    uCellI_ = nCellsI_ - 1;
    uCellJ_ = nCellsJ_ - 1;
    uCellK_ = nCellsK_ - 1;

    uFaceI_ = nFacesI_ - 1;
    uFaceJ_ = nFacesJ_ - 1;
    uFaceK_ = nFacesK_ - 1;

    conservedFieldName_ = conservedFieldName;
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
                    gradPhiField(i, j, k) = (phiField.faceE(i, j, k)*mesh.fAreaNormE(i, j, k)
                                             + phiField.faceW(i, j, k)*mesh.fAreaNormW(i, j, k)
                                             + phiField.faceN(i, j, k)*mesh.fAreaNormN(i, j, k)
                                             + phiField.faceS(i, j, k)*mesh.fAreaNormS(i, j, k)
                                             + phiField.faceT(i, j, k)*mesh.fAreaNormT(i, j, k)
                                             + phiField.faceB(i, j, k)*mesh.fAreaNormB(i, j, k))/mesh.cellVol(i, j, k);
                }
            }
        }

        break;
    };
}

void FvScheme::computeCellCenteredGradients(Field<Vector3D> &vecField, Field<Tensor3D> &gradVecField, int method)
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
                        gradVecField(i, j, k)(l, m) = x(m, 0);
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
                    gradVecField(i, j, k) = (tensor(vecField.faceE(i, j, k), mesh.fAreaNormE(i, j, k))
                                            + tensor(vecField.faceW(i, j, k), mesh.fAreaNormW(i, j, k))
                                            + tensor(vecField.faceN(i, j, k), mesh.fAreaNormN(i, j, k))
                                            + tensor(vecField.faceS(i, j, k), mesh.fAreaNormS(i, j, k))
                                            + tensor(vecField.faceT(i, j, k), mesh.fAreaNormT(i, j, k))
                                            + tensor(vecField.faceB(i, j, k), mesh.fAreaNormB(i, j, k)))/mesh.cellVol(i, j, k);
                }
            }
        }

        break;
    };
}

void FvScheme::computeFaceCenteredGradients(Field<double> &phiField, Field<Vector3D> &gradPhiField)
{

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
