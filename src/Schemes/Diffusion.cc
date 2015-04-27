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

void Diffusion::computeFaceFluxes()
{
    int i, j, k, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& phiField = *phiFieldPtr_;

    uI = nFacesI_ - 1;
    uJ = nFacesJ_ - 1;
    uK = nFacesK_ - 1;

    for(k = 0; k < nFacesK_; ++k)
    {
        for(j = 0; j < nFacesJ_; ++j)
        {
            for(i = 0; i < nFacesI_; ++i)
            {
                if(j < uJ && k < uK)
                    phiField.fluxI(i, j, k) = dot(gradPhiField_.faceI(i, j, k), mesh.fAreaNormI(i, j, k));

                if(i < uI && k < uK)
                    phiField.fluxJ(i, j, k) = dot(gradPhiField_.faceJ(i, j, k), mesh.fAreaNormJ(i, j, k));

                if(i < uI && j < uJ)
                    phiField.fluxK(i, j, k) = dot(gradPhiField_.faceK(i, j, k), mesh.fAreaNormK(i, j, k));
            }
        }
    }
}

// ************* Public Methods *************

void Diffusion::initialize(Input &input, HexaFvmMesh &mesh, std::string conservedFieldName)
{
    FvScheme::initialize(input, mesh, conservedFieldName);
    phiFieldPtr_ = &mesh.findScalarField(conservedFieldName_);
    mesh.addVectorField("phiGrad");

    gradPhiField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
}

int Diffusion::nConservedVariables()
{
    return phiFieldPtr_->size();
}

void Diffusion::discretize(std::vector<double>& timeDerivatives)
{
    int i, j, k, l;
    Field<double>& phiField = *phiFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;

    phiField.setBoundaryFields();
    computeCellCenteredGradients(phiField, gradPhiField_);
    computeFaceCenteredGradients(phiField, gradPhiField_);
    computeFaceFluxes();

    for(k = 0, l = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
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
