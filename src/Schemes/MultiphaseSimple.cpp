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

#include <iomanip>

#include "MultiphaseSimple.h"

MultiphaseSimple::MultiphaseSimple()
    :
      alphaField0_("alphaField0", AUXILLARY),
      gradAlphaField_("gradAlphaField", PRIMITIVE),
      interfaceNormals_("interfaceNormals", AUXILLARY),
      kField_("k", AUXILLARY),
      bF_("bF", AUXILLARY),
      sigma_(0.073)
{

}

void MultiphaseSimple::computePhysicalConstants(Field<double> &alphaField, Field<double> &rhoField, Field<double>& muField)
{
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                rhoField(i, j, k) = (1. - alphaField(i, j, k))*rho1_ + alphaField(i, j, k)*rho2_;
                muField(i, j, k) = (1. - alphaField(i, j, k))*mu1_ + alphaField(i, j, k)*mu2_;
            }
        }
    }
}

void MultiphaseSimple::computeCurvature(Field<double> &alphaField)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;
    double mag;

    extrapolateInteriorFaces(alphaField, gradAlphaField_);
    computeCellCenteredGradients(alphaField, gradAlphaField_, DIVERGENCE_THEOREM);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                mag = gradAlphaField_(i, j, k).mag();

                if(mag < 1e-15)
                    interfaceNormals_(i, j, k) = Vector3D(0., 0., 0.);
                else
                    interfaceNormals_(i, j, k) = gradAlphaField_(i, j, k)/mag;
            }
        }
    }

    extrapolateInteriorFaces(interfaceNormals_, gradVecField_);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                kField_(i, j, k) = -(dot(interfaceNormals_.faceE(i, j, k), mesh.fAreaNormE(i, j, k))
                                     + dot(interfaceNormals_.faceW(i, j, k), mesh.fAreaNormW(i, j, k))
                                     + dot(interfaceNormals_.faceN(i, j, k), mesh.fAreaNormN(i, j, k))
                                     + dot(interfaceNormals_.faceS(i, j, k), mesh.fAreaNormS(i, j, k))
                                     + dot(interfaceNormals_.faceT(i, j, k), mesh.fAreaNormT(i, j, k))
                                     + dot(interfaceNormals_.faceB(i, j, k), mesh.fAreaNormB(i, j, k)))/mesh.cellVol(i, j, k);
            }
        }
    }
}

void MultiphaseSimple::computeSurfaceTensionSource()
{
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                bF_(i, j, k) = sigma_*kField_(i, j, k)*gradAlphaField_(i, j, k);
            }
        }
    }
}

void MultiphaseSimple::advectAlphaField(Field<double>& rhoField, Field<double>& massFlowField, Field<Vector3D> &uField, double timeStep, Field<double> &alphaField)
{

}

void MultiphaseSimple::initialize(Input &input, HexaFvmMesh& mesh)
{
    Simple::initialize(input, mesh);

    alphaFieldPtr_ = &mesh.findScalarField("alpha");

    alphaField0_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradAlphaField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    interfaceNormals_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    kField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    bF_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    rho1_ = input.inputDoubles["rho"];
    mu1_ = input.inputDoubles["mu"];
    rho2_ = input.inputDoubles["rho2"];
    mu2_ = input.inputDoubles["mu2"];
}

void MultiphaseSimple::storeAlphaField(Field<double> &alphaField)
{
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                alphaField0_(i, j, k) = alphaField(i, j, k);
            }
        }
    }
}

void MultiphaseSimple::discretize(double timeStep, std::vector<double> &timeDerivatives)
{
    Field<double>& rhoField = *rhoFieldPtr_;
    Field<double>& muField = *muFieldPtr_;
    Field<double>& massFlowField = *massFlowFieldPtr_;
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    Field<double>& alphaField = *alphaFieldPtr_;
    int i;

    storeUField(uField, uField0_);
    //storeAlphaField(alphaField);

    //computePhysicalConstants(alphaField, rhoField, muField);

    for(i = 0; i < maxInnerIters_; ++i)
    {
        //computeCurvature(alphaField);
        computeMomentum(rhoField, muField, massFlowField, &bF_, timeStep, uField, pField);
        computePCorr(rhoField, massFlowField, uField, pField);
        correctContinuity(rhoField, massFlowField, uField, pField);
        std::cout << "\rDTS iteration completion  |      " << (i + 1) << "/" << maxInnerIters_ << std::fixed << std::setprecision(2) << " (" << 100.*(i + 1)/maxInnerIters_ << "%)";
    }

    //advectAlphaField(rhoField, massFlowField, uField, timeStep, alphaField);
}
