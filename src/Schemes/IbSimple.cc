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

    //- Determine the solid cells, set the rest as fluid
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(ibSphere_.isInside(mesh.cellXc(i, j, k)))
                {
                    cellStatus_(i, j, k) = INACTIVE;
                    ibField_(i, j, k) = SOLID;
                    uField(i, j, k) = Vector3D(0., 0., 0.);
                }
                else
                {
                    cellStatus_(i, j, k) = ACTIVE;
                    ibField_(i, j, k) = FLUID;
                }

                mesh.findScalarField("ibField")(i, j, k) = ibField_(i, j, k);
            }
        }
    }

    //- Determine the immersed boundary cells
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(ibField_(i, j, k) == SOLID)
                    continue;

                if(ibField_(i + 1, j, k) == SOLID ||
                        ibField_(i - 1, j, k) == SOLID ||
                        ibField_(i, j + 1, k) == SOLID ||
                        ibField_(i, j - 1, k) == SOLID ||
                        ibField_(i, j, k + 1) == SOLID ||
                        ibField_(i, j, k - 1) == SOLID)
                {
                    ibField_(i, j, k) = IB;
                    cellStatus_(i, j, k) = INTERPOLATION;
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
    int i, j, k, l, n;
    Point3D tmpPoints[6];
    double tmpValues[6];

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(ibField_(i, j, k) == IB)
                {
                    tmpPoints[0] = mesh.cellXc(i + 1, j, k);
                    tmpPoints[1] = mesh.cellXc(i - 1, j, k);
                    tmpPoints[2] = mesh.cellXc(i, j + 1, k);
                    tmpPoints[3] = mesh.cellXc(i, j - 1, k);
                    tmpPoints[4] = mesh.cellXc(i, j, k + 1);
                    tmpPoints[5] = mesh.cellXc(i, j, k - 1);
                    for(l = 0; l < 3; ++l)
                    {
                        tmpValues[0] = uField(i + 1, j, k)(l);
                        tmpValues[1] = uField(i - 1, j, k)(l);
                        tmpValues[2] = uField(i, j + 1, k)(l);
                        tmpValues[3] = uField(i, j - 1, k)(l);
                        tmpValues[4] = uField(i, j, k + 1)(l);
                        tmpValues[5] = uField(i, j, k - 1)(l);

                        uField(i, j, k)(l) = Interpolation::linear(tmpPoints, tmpValues, 6, mesh.cellXc(i, j, k));
                    }

                    n = 0;
                    pField(i, j, k) = 0.;
                    pCorr_(i, j, k) = 0.;

                    if(ibField_(i + 1, j, k) == FLUID)
                    {
                        pField(i, j, k) += pField(i + 1, j, k);// + dot(gradPField_(i + 1, j, k), mesh.rCellW(i + 1, j, k));
                        pCorr_(i, j, k) += pCorr_(i + 1, j, k);// + dot(gradPCorr_(i + 1, j, k), mesh.rCellW(i + 1, j, k));
                        ++n;
                    }

                    if(ibField_(i - 1, j, k) == FLUID)
                    {
                        pField(i, j, k) += pField(i - 1, j, k);// + dot(gradPField_(i - 1, j, k), mesh.rCellE(i - 1, j, k));
                        pCorr_(i, j, k) += pCorr_(i - 1, j, k);// + dot(gradPCorr_(i - 1, j, k), mesh.rCellE(i - 1, j, k));
                        ++n;
                    }

                    if(ibField_(i, j + 1, k) == FLUID)
                    {
                        pField(i, j, k) += pField(i, j + 1, k);// + dot(gradPField_(i, j + 1, k), mesh.rCellS(i, j + 1, k));
                        pCorr_(i, j, k) += pCorr_(i, j + 1, k);// + dot(gradPCorr_(i, j + 1, k), mesh.rCellS(i, j + 1, k));
                        ++n;
                    }

                    if(ibField_(i, j - 1, k) == FLUID)
                    {
                        pField(i, j, k) += pField(i, j - 1, k);/// + dot(gradPField_(i, j - 1, k), mesh.rCellN(i, j - 1, k));
                        pCorr_(i, j, k) += pCorr_(i, j - 1, k);// + dot(gradPCorr_(i, j - 1, k), mesh.rCellN(i, j - 1, k));
                        ++n;
                    }

                    if(ibField_(i, j, k + 1) == FLUID)
                    {
                        pField(i, j, k) += pField(i, j, k + 1);// + dot(gradPField_(i, j, k + 1), mesh.rCellB(i, j, k + 1));
                        pCorr_(i, j, k) += pCorr_(i, j, k + 1);// + dot(gradPCorr_(i, j, k + 1), mesh.rCellB(i, j, k + 1));
                        ++n;
                    }

                    if(ibField_(i, j, k - 1) == FLUID)
                    {
                        pField(i, j, k) += pField(i, j, k - 1);// + dot(gradPField_(i, j, k - 1), mesh.rCellT(i, j, k - 1));
                        pCorr_(i, j, k) += pCorr_(i, j, k - 1);// + dot(gradPCorr_(i, j, k - 1), mesh.rCellT(i, j, k - 1));
                        ++n;
                    }

                    pField(i, j, k) /= double(n);
                    pCorr_(i, j, k) /= double(n);
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
    }
}
