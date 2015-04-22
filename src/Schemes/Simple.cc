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
      gradUField_("gradUField", PRIMITIVE),
      gradPField_("gradPField", PRIMITIVE),
      relaxationFactor_(0.8),
      rho_(998.),
      mu_(0.1)
{
    nu_ = mu_/rho_;
}

// ************* Private Methods *************

// ************* Public Methods *************

void Simple::initialize(Input &input, HexaFvmMesh &mesh)
{
    int nCellsI, nCellsJ, nCellsK;

    FvScheme::initialize(input, mesh, "NA");
    uFieldPtr_ = &mesh.findVectorField("u");
    pFieldPtr_ = &mesh.findScalarField("p");

    relaxationFactor_ = input.inputDoubles["relaxationFactor"];
    rho_ = input.inputDoubles["rho"];
    mu_ = input.inputDoubles["mu"];
    nu_ = mu_/rho_;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    gradUField_.allocate(nCellsI, nCellsJ, nCellsK);
    gradPField_.allocate(nCellsI, nCellsJ, nCellsK);
}

int Simple::nConservedVariables()
{
    return 4*meshPtr_->size();
}

void Simple::discretize(std::vector<double> &timeDerivatives_)
{
    int i, j, k, l, nCellsI, nCellsJ, nCellsK;
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    //- Set the boundary fields

    uField.setBoundaryFields();
    pField.setBoundaryFields();

    //- Compute relevant gradients and jacobians

    computeCellCenteredJacobians(uField, gradUField_);
    computeCellCenteredGradients(pField, gradPField_);

    computeUpwindFaceCenteredReconstruction(uField, gradUField_, uField);
    computeFaceCenteredGradients(pField, gradPField_);

    //- Compute momentum using currently available pressure field

    computeMomentum(uField, gradUField_, gradPField_);

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

}

void Simple::updateSolution(std::vector<double> &update, int method)
{

}

void Simple::computeMomentum(Field<Vector3D> &uField, Field<Tensor3D> &gradUField, Field<Vector3D> &gradPField)
{
    int i, j, k, nCellsI, nCellsJ, nCellsK;
    HexaFvmMesh& mesh = *meshPtr_;

    nCellsI = mesh.nCellsI();
    nCellsJ = mesh.nCellsJ();
    nCellsK = mesh.nCellsK();

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                //- Convective flux

                uField.fluxE(i, j, k) = -(tensor(uField.faceE(i, j, k), uField.faceE(i, j, k))*mesh.fAreaNormE(i, j, k));
                uField.fluxN(i, j, k) = -(tensor(uField.faceN(i, j, k), uField.faceN(i, j, k))*mesh.fAreaNormN(i, j, k));
                uField.fluxT(i, j, k) = -(tensor(uField.faceT(i, j, k), uField.faceT(i, j, k))*mesh.fAreaNormT(i, j, k));

                //- Diffusive flux

                uField.fluxE(i, j, k) += nu_*gradUField.faceE(i, j, k)*mesh.fAreaNormE(i, j, k);
                uField.fluxN(i, j, k) += nu_*gradUField.faceN(i, j, k)*mesh.fAreaNormN(i, j, k);
                uField.fluxT(i, j, k) += nu_*gradUField.faceT(i, j, k)*mesh.fAreaNormT(i, j, k);
            }
        }
    }
}
