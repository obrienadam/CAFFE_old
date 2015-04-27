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

// ************* Constructors and Destructors *************

Simple::Simple()
    :
      massFlow_("massFlow", PRIMITIVE),
      pCorr_("pCorr", PRIMITIVE),
      gradUField_("gradUField", PRIMITIVE),
      gradPField_("gradPField", PRIMITIVE),
      relaxationFactor_(0.8),
      rho_(998.),
      mu_(0.1)
{
    nu_ = mu_/rho_;
}

// ************* Private Methods *************

void Simple::computeMassFlux()
{
    int i, j, k, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    double alpha;
    Vector3D uBar, gradP, gradPBar;
    Tensor3D dBar;

    //- Set the east and west mass fluxes at the boundaries

    uI = nFacesI_ - 1;

    for(k = 0; k < nFacesK_; ++k)
    {
        for(j = 0; j < nFacesJ_; ++j)
        {
            massFlow_.faceI(uI, j, k) = dot(rho_*uField.faceI(uI, j, k), mesh.fAreaNormI(uI, j, k));
            massFlow_.faceI(0, j, k) = dot(rho_*uField.faceI(0, j, k), mesh.fAreaNormI(0, j, k));
        }
    }

    //- Set the north and south mass fluxes at the boundaries

    uJ = nFacesJ_ - 1;

    for(k = 0; k < nFacesK_; ++k)
    {
        for(i = 0; i < nFacesI_; ++i)
        {
            massFlow_.faceJ(i, uJ, k) = dot(rho_*uField.faceJ(i, uJ, k), mesh.fAreaNormJ(i, uJ, k));
            massFlow_.faceJ(i, 0, k) = dot(rho_*uField.faceJ(i, 0, k), mesh.fAreaNormJ(i, 0, k));
        }
    }

    //- Set the top and bottom mass fluxes at the boundaries

    uK = nFacesK_ - 1;

    for(j = 0; j < nFacesJ_; ++j)
    {
        for(i = 0; i < nFacesI_; ++i)
        {
            massFlow_.faceK(i, j, uK) = dot(rho_*uField.faceK(i, j, uK), mesh.fAreaNormK(i, j, uK));
            massFlow_.faceK(i, j, 0) = dot(rho_*uField.faceK(i, j, 0), mesh.fAreaNormK(i, j, 0));
        }
    }

    uI = nCellsI_ - 1;
    uJ = nCellsJ_ - 1;
    uK = nCellsK_ - 1;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(i < uI)
                {
                    alpha = getAlpha(i, j, k, EAST);

                    uBar = alpha*uField(i, j, k) + (1. - alpha)*uField(i + 1, j, k);
                    gradPBar = alpha*gradPField_(i, j, k) + (1. - alpha)*gradPField_(i + 1, j, k);
                    dBar = alpha*dField_(i, j, k) + (1. - alpha)*dField_(i + 1, j, k);
                    gradP = (pField(i + 1, j, k) - pField(i, j, k))/mesh.rCellMagE(i, j, k)*mesh.rnCellE(i, j, k);

                    massFlow_.faceE(i, j, k) = dot(rho_*uBar - rho_*dBar*(gradP - gradPBar), mesh.fAreaNormE(i, j, k));
                }

                if(j < uJ)
                {
                    alpha = getAlpha(i, j, k, NORTH);

                    uBar = alpha*uField(i, j, k) + (1. - alpha)*uField(i, j + 1, k);
                    gradPBar = alpha*gradPField_(i, j, k) + (1. - alpha)*gradPField_(i, j + 1, k);
                    dBar = alpha*dField_(i, j, k) + (1. - alpha)*dField_(i, j + 1, k);
                    gradP = (pField(i, j + 1, k) - pField(i, j, k))/mesh.rCellMagN(i, j, k)*mesh.rnCellN(i, j, k);

                    massFlow_.faceN(i, j, k) = dot(rho_*uBar - rho_*dBar*(gradP - gradPBar), mesh.fAreaNormN(i, j, k));
                }

                if(k < uK)
                {
                    alpha = getAlpha(i, j, k, TOP);

                    uBar = alpha*uField(i, j, k) + (1. - alpha)*uField(i, j, k + 1);
                    gradPBar = alpha*gradPField_(i, j, k) + (1. - alpha)*gradPField_(i, j, k + 1);
                    dBar = alpha*dField_(i, j, k) + (1. - alpha)*dField_(i, j, k + 1);
                    gradP = (pField(i, j, k + 1) - pField(i, j, k))/mesh.rCellMagT(i, j, k)*mesh.rnCellT(i, j, k);

                    massFlow_.faceT(i, j, k) = dot(rho_*uBar - rho_*dBar*(gradP - gradPBar), mesh.fAreaNormT(i, j, k));
                }
            }
        }
    }
}

// ************* Public Methods *************

void Simple::initialize(Input &input, HexaFvmMesh &mesh)
{
    FvScheme::initialize(input, mesh, "NA");
    uFieldPtr_ = &mesh.findVectorField("u");
    pFieldPtr_ = &mesh.findScalarField("p");

    relaxationFactor_ = input.inputDoubles["relaxationFactor"];
    rho_ = input.inputDoubles["rho"];
    mu_ = input.inputDoubles["mu"];
    nu_ = mu_/rho_;

    massFlow_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    pCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradUField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradPField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
}

int Simple::nConservedVariables()
{
    return 4*meshPtr_->size();
}

void Simple::discretize(std::vector<double> &timeDerivatives_)
{
    int i, j, k, l;
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;

    //- Set the boundary fields

    uField.setBoundaryFields();
    pField.setBoundaryFields();

    //- Compute relevant gradients and jacobians

    computeCellCenteredJacobians(uField, gradUField_);
    computeCellCenteredGradients(pField, gradPField_);

    //- Use the Rhie-Chow interpolation to compute the mass fluxes

    computeMassFlux();

    //- Compute the predicted momentum
}

void Simple::copySolution(std::vector<double> &original)
{

}

void Simple::updateSolution(std::vector<double> &update, int method)
{

}
