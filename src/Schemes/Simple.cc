/**
 * @file    Simple.cc
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
 * Simple.
 */

#include "Simple.h"

void Simple::initialize(Input &input, HexaFvmMesh &mesh)
{
    int nCellsI, nCellsJ, nCellsK;

    FvScheme::initialize(input, mesh, "NA");
    uFieldPtr_ = &mesh.findScalarField("u");
    vFieldPtr_ = &mesh.findScalarField("v");
    wFieldPtr_ = &mesh.findScalarField("w");
    pFieldPtr_ = &mesh.findScalarField("p");

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    gradUField_.allocate(nCellsI, nCellsJ, nCellsK);
    gradVField_.allocate(nCellsI, nCellsJ, nCellsK);
    gradWField_.allocate(nCellsI, nCellsJ, nCellsK);
    gradPField_.allocate(nCellsI, nCellsJ, nCellsK);
}

int Simple::nConservedVariables()
{
    return 4*meshPtr_->size();
}

void Simple::discretize(std::vector<double> &timeDerivatives_)
{
    int i, j, k, l, nCellsI, nCellsJ, nCellsK;
    Field<double>& uField = *uFieldPtr_;
    Field<double>& vField = *vFieldPtr_;
    Field<double>& wField = *wFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    //- Set the boundary fields

    uField.setBoundaryFields();
    vField.setBoundaryFields();
    wField.setBoundaryFields();
    pField.setBoundaryFields();

    computeCellCenteredGradients(uField, gradUField_);
    computeCellCenteredGradients(vField, gradVField_);
    computeCellCenteredGradients(wField, gradWField_);
    computeCellCenteredGradients(pField, gradPField_);

    for(k = 0, l = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                //timeDerivatives[l] = phiField.sumFluxes(i, j, k)/mesh.cellVol(i, j, k);
                //mesh.findVectorField("phiGrad")(i, j, k) = gradPhiField_(i, j, k);
                //++l;
            }
        }
    }
}

void Simple::copySolution(std::vector<double> &original)
{
    int k, l, size;
    Field<double>& uField = *uFieldPtr_;
    Field<double>& vField = *vFieldPtr_;
    Field<double>& wField = *wFieldPtr_;
    Field<double>& pField = *pFieldPtr_;

    size = uField.size();

    for(k = 0, l = 0; k < size; ++k, ++l)
    {
        original[l] = uField(k);
    }

    for(k = 0; k < size; ++k, ++l)
    {
        original[l] = vField(k);
    }

    for(k = 0; k < size; ++k, ++l)
    {
        original[l] = wField(k);
    }

    for(k = 0; k < size; ++k, ++l)
    {
        original[l] = pField(k);
    }
}

void Simple::updateSolution(std::vector<double> &update, int method)
{
    int k, l, size;
    Field<double>& uField = *uFieldPtr_;
    Field<double>& vField = *vFieldPtr_;
    Field<double>& wField = *wFieldPtr_;
    Field<double>& pField = *pFieldPtr_;

    size = uField.size();

    switch (method)
    {
    case ADD:

        for(k = 0, l = 0; k < size; ++k, ++l)
        {
            uField(k) += update[l];
        }

        for(k = 0; k < size; ++k, ++l)
        {
            vField(k) += update[l];
        }

        for(k = 0; k < size; ++k, ++l)
        {
            wField(k) += update[l];
        }

        for(k = 0; k < size; ++k, ++l)
        {
            pField(k) += update[l];
        }

        break;
    case REPLACE:

        for(k = 0, l = 0; k < size; ++k, ++l)
        {
            uField(k) = update[l];
        }

        for(k = 0; k < size; ++k, ++l)
        {
            vField(k) = update[l];
        }

        for(k = 0; k < size; ++k, ++l)
        {
            wField(k) = update[l];
        }

        for(k = 0; k < size; ++k, ++l)
        {
            pField(k) = update[l];
        }

        break;
    };
}
