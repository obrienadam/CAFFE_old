/**
 * @file    IbPiso.cpp
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

#include "IbPiso.h"
#include "Interpolation.h"

/***************** Constructors and destructors *****************************/

IbPiso::IbPiso(const Input &input, const HexaFvmMesh &mesh)
    :
      Piso(input, mesh),
      ibField_(mesh, Field<CellType>::AUXILLARY, "ibField")
{
    ibSphere_.radius = input.caseParameters.get<double>("ImmersedBoundaries.Sphere.radius");
    ibSphere_.center = std::stov(input.caseParameters.get<std::string>("ImmersedBoundaries.Sphere.center"));
}

/***************** Public methods *****************************/

double IbPiso::solve(double timeStep)
{
    int i, j;

    uField0_ = uField_;
    computeIbField();

    for(j = 0; j < nInnerIters_; ++j)
    {
        computeIbCoeffs();
        computeMomentum(timeStep);

        for(i = 0; i < nPCorrections_; ++i)
        {
            computePCorr();
            correct();
        }
    }

    return computeContinuityError();
}

/***************** Private methods *****************************/

void IbPiso::computeIbField()
{
    int i, j, k;

    //- This procedure is done as per Mittal et al. (2008)
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(ibSphere_.isInside(mesh_.cellXc(i, j, k)))
                {
                    ibField_(i, j, k) = SOLID;
                    mesh_.iMap.setInactive(i, j, k);
                    uField_(i, j, k) = Vector3D(0., 0., 0.);
                    pField_(i, j, k) = 0.;
                }
                else
                {
                    ibField_(i, j, k) = FLUID;
                    mesh_.iMap.setActive(i, j, k);
                }
            }
        }
    }

    //- Determine the IB or "ghost" cells
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(ibField_(i, j, k) == SOLID)
                {
                    if(ibField_(i + 1, j, k) == FLUID || ibField_(i - 1, j, k) == FLUID
                            || ibField_(i, j + 1, k) == FLUID || ibField_(i, j - 1, k) == FLUID
                            || ibField_(i, j, k + 1) == FLUID || ibField_(i, j, k - 1) == FLUID)
                    {
                        ibField_(i, j, k) = IB;
                        mesh_.iMap.setGhost(i, j, k);
                    }
                }
            }
        }
    }
    mesh_.iMap.generateGlobalIndices();

    destroyMatrices();
    createMatrices(2, 3, 9);
}

void IbPiso::computeIbCoeffs()
{
    using namespace std;

    Point3D boundaryPoint, imagePoint, tmpPoints[8];
    int i, j, k, m, ii[8], jj[8], kk[8], pointNo, colNos[9];
    Matrix beta(1, 8);
    double values[9];

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(ibField_(i, j, k) == IB)
                {
                    // Find the boundary and image points
                    boundaryPoint = ibSphere_.nearestIntersect(mesh_.cellXc(i, j, k));
                    imagePoint = 2.*boundaryPoint - mesh_.cellXc(i, j , k);

                    mesh_.locateEnclosingCells(imagePoint, ii, jj, kk);

                    for(pointNo = 0; pointNo < 8; ++pointNo)
                        tmpPoints[pointNo] = mesh_.cellXc(ii[pointNo], jj[pointNo], kk[pointNo]);

                    //- Determine the interpolation coefficients for the image point
                    beta = Interpolation::computeTrilinearCoeffs(tmpPoints, 8, imagePoint);

                    values[0] = 1.;
                    colNos[0] = mesh_.iMap(i, j, k, 0);

                    for(m = 0; m < 8; ++m)
                    {
                        if(ii[m] == i && jj[m] == j && kk[m] == k)
                        {
                            values[0] += beta(0, m);
                            colNos[m + 1] = -1;
                        }
                        else
                        {
                            values[m + 1] = beta(0, m);
                            colNos[m + 1] = mesh_.iMap(ii[m], jj[m], kk[m], 0);
                        }
                    }
                    A_[0].setRow(mesh_.iMap(i, j, k, 0), 9, colNos, values);

                    //- Add the equation for the ghost cells to the coefficient matrix form pressure
                    for(m = 0; m < 8; ++m)
                    {
                        values[0] = -1;
                        colNos[0] = mesh_.iMap(i, j, k, 0);

                        if(ii[m] == i && jj[m] == j && kk[m] == k)
                        {
                            values[0] += beta(0, m);
                            colNos[m + 1] = -1;
                        }
                        else
                        {
                            values[m + 1] = beta(0, m);
                            colNos[m + 1] = mesh_.iMap(ii[m], jj[m], kk[m], 0);
                        }
                    }

                    A_[1].setRow(mesh_.iMap(i, j, k, 0), 9, colNos, values);
                }
            }
        }
    }
}
