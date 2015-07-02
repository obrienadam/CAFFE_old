/**
 * @file    IbSimple.cc
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
 * IbSimple.
 */

#include <iomanip>

#include "IbSimple.h"
#include "Interpolation.h"

IbSimple::IbSimple()
    :
      ibField_("ib", PRIMITIVE),
      ibSourceField_("ibSource", AUXILLARY)
{

}

void IbSimple::computeIbField(Field<Vector3D>& uField, Field<double>& pField)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;

    //- This procedure is done as per Mittal et al. (2008)
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(ibSphere_.isInside(mesh.cellXc(i, j, k)))
                {
                    ibField_(i, j, k) = SOLID;
                    cellStatus_(i, j, k) = INACTIVE;
                }
                else
                {
                    ibField_(i, j, k) = FLUID;
                    cellStatus_(i, j, k) = ACTIVE;
                }
            }
        }
    }

    //- Determine the IB or "ghost" cells
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(ibField_(i, j, k) == SOLID)
                {
                    if(ibField_(i + 1, j, k) == FLUID || ibField_(i - 1, j, k) == FLUID
                            || ibField_(i, j + 1, k) == FLUID || ibField_(i, j - 1, k) == FLUID
                            || ibField_(i, j, k + 1) == FLUID || ibField_(i, j, k - 1) == FLUID)
                    {
                        ibField_(i, j, k) = IB;
                        cellStatus_(i, j, k) = GHOST;
                    }
                }

                mesh.findScalarField("ibField")(i, j, k) = ibField_(i, j, k);
            }
        }
    }
}

void IbSimple::setIbCells(Field<Vector3D>& uField, Field<double>& pField)
{
    using namespace std;

    HexaFvmMesh& mesh = *meshPtr_;
    Point3D boundaryPoint, imagePoint, tmpPoints[8];
    int i, j, k, ii[8], jj[8], kk[8], pointNo;
    Matrix beta(1, 8);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(ibField_(i, j, k) == IB)
                {
                    // Find the boundary and image points
                    boundaryPoint = ibSphere_.nearestIntersect(mesh.cellXc(i, j, k));
                    imagePoint = 2.*boundaryPoint - mesh.cellXc(i, j , k);

                    mesh.locateEnclosingCells(imagePoint, ii, jj, kk);

                    for(pointNo = 0; pointNo < 8; ++pointNo)
                        tmpPoints[pointNo] = mesh.cellXc(ii[pointNo], jj[pointNo], kk[pointNo]);

                    beta = Interpolation::computeTrilinearCoeffs(tmpPoints, 8, imagePoint);
                }
            }
        }
    }
}

void IbSimple::initialize(Input &input, HexaFvmMesh &mesh)
{
    Simple::initialize(input, mesh);

    ibSphere_.radius = input.inputDoubles["ibSphereRadius"];
    ibSphere_.center = stov(input.inputStrings["ibSphereCenter"]);

    ibField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    ibSourceField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
}

void IbSimple::discretize(double timeStep, std::vector<double> &timeDerivatives)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    Field<double>& rhoField = *rhoFieldPtr_;
    Field<double>& muField = *muFieldPtr_;
    Field<double>& massFlowField = *massFlowFieldPtr_;
    int i;

    storeUField(uField, uField0_);
    computeIbField(uField, pField);

    for(i = 0; i < maxInnerIters_; ++i)
    {
        setIbCells(uField, pField);
        computeMomentum(rhoField, muField, massFlowField, NULL, timeStep, uField, pField);
        computePCorr(rhoField, massFlowField, uField, pField);
        correctContinuity(rhoField, massFlowField, uField, pField);
        std::cout << "\rDTS iteration completion  |      " << (i + 1) << "/" << maxInnerIters_ << std::fixed << std::setprecision(2) << " (" << 100.*(i + 1)/maxInnerIters_ << "%)";
    }
}
