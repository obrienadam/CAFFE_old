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

IbSimple::IbSimple()
    :
      ibField_("ib", AUXILLARY),
      ibSourceField_("ibSource", AUXILLARY)
{

}

bool IbSimple::isSolidCell(int i, int j, int k, HexaFvmMesh& mesh)
{
    if(ibSphere_.isInside(mesh.cellXc(i, j, k)))
        return true;
    else
        return false;
}

bool IbSimple::isFluidCell(int i, int j, int k, HexaFvmMesh &mesh)
{
    for(int l = 0; i < 8; ++i)
    {
        if(ibSphere_.isInside(mesh.node(i, j, k, l)))
            return false;
    }

    return true;
}

void IbSimple::computeIbField()
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(isFluidCell(i, j, k, mesh))
                {
                    ibField_(i, j, k) = FLUID;
                    cellStatus_(i, j, k) = ACTIVE;
                }
                else if (isSolidCell(i, j, k, mesh))
                {
                    ibField_(i, j, k) = SOLID;
                    cellStatus_(i, j, k) = INACTIVE;
                }
                else
                {
                    ibField_(i, j, k) = IB;
                    cellStatus_(i, j, k) = BOUNDARY;
                }
            }
        }
    }
}

void IbSimple::computeIbSourceTerm()
{

}

void IbSimple::initialize(Input &input, HexaFvmMesh &mesh)
{
    Simple::initialize(input, mesh);

    ibSphere_.radius = input.inputDoubles["ibSphereRadius"];
    ibSphere_.center = Point3D(input.inputStrings["ibSphereCenter"]);

    ibField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    ibSourceField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
}

void IbSimple::discretize(double timeStep, std::vector<double> &timeDerivatives)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    Field<double>& rhoField = *rhoFieldPtr_;
    Field<double>& muField = *muFieldPtr_;
    int i;

    storeUField(uField, uField0_);

    for(i = 0; i < maxInnerItrs_; ++i)
    {
        computeIbField();
        computeIbSourceTerm();
        computeMomentum(rhoField, muField, NULL, timeStep, uField, pField);
        computePCorr(rhoField, uField, pField);
        correctContinuity(rhoField, uField, pField);
    }
}
