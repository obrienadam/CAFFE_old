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
    using namespace std;

    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;
    SparseMatrix A;
    SparseVector x, b;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                if(timeAccurate_)
                    a0P_(i, j, k) = mesh.cellVol(i, j, k)/timeStep;
                else
                    a0P_(i, j, k) = 0.;

                // Face coefficients
                aE_(i, j, k) =  min(massFlowField.faceE(i, j, k), 0.)/rhoField.faceE(i, j, k);
                aW_(i, j, k) =  -max(massFlowField.faceW(i, j, k), 0.)/rhoField.faceW(i, j, k);
                aN_(i, j, k) =  min(massFlowField.faceN(i, j, k), 0.)/rhoField.faceN(i, j, k);
                aS_(i, j, k) =  -max(massFlowField.faceS(i, j, k), 0.)/rhoField.faceS(i, j, k);
                aT_(i, j, k) =  min(massFlowField.faceT(i, j, k), 0.)/rhoField.faceT(i, j, k);
                aB_(i, j, k) =  -max(massFlowField.faceB(i, j, k), 0.)/rhoField.faceB(i, j, k);

                // Central coefficient
                aP_(i, j, k) = max(massFlowField.faceE(i, j, k), 0.)/rhoField.faceE(i, j, k)
                        - min(massFlowField.faceW(i, j, k), 0.)/rhoField.faceW(i, j, k)
                        + max(massFlowField.faceN(i, j, k), 0.)/rhoField.faceN(i, j, k)
                        - min(massFlowField.faceS(i, j, k), 0.)/rhoField.faceS(i, j, k)
                        + max(massFlowField.faceT(i, j, k), 0.)/rhoField.faceT(i, j, k)
                        - min(massFlowField.faceB(i, j, k), 0.)/rhoField.faceB(i, j, k)
                        + a0P_(i, j, k);

                bP_(i, j, k).x = a0P_(i, j, k)*alphaField0_(i, j, k);

                //- Boundary conditions
                // I-direction bcs
                if(i == uCellI_ && uField.getEastBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aE_(i, j, k);
                    aE_(i, j, k) = 0.;
                }
                else if(i == uCellI_ || cellStatus_(i + 1, j, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aE_(i, j, k)*alphaField(i + 1, j, k);
                    aE_(i, j, k) = 0.;
                }

                if(i == 0 && uField.getWestBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aW_(i, j, k);
                    aW_(i, j, k) = 0.;
                }
                else if(i == 0 || cellStatus_(i - 1, j, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aW_(i, j, k)*alphaField(i - 1, j, k);
                    aW_(i, j, k) = 0.;
                }

                // J-direction bcs
                if(j == uCellJ_ && uField.getNorthBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aN_(i, j, k);
                    aN_(i, j, k) = 0.;
                }
                else if(j == uCellJ_ || cellStatus_(i, j + 1, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aN_(i, j, k)*alphaField(i, j + 1, k);
                    aN_(i, j, k) = 0.;
                }

                if(j == 0 && uField.getSouthBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aS_(i, j, k);
                    aS_(i, j, k) = 0.;
                }
                else if(j == 0 || cellStatus_(i, j - 1, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aS_(i, j, k)*alphaField(i, j - 1, k);
                    aS_(i, j, k) = 0.;
                }

                // K-direction bcs
                if(k == uCellK_ && uField.getTopBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aT_(i, j, k);
                    aT_(i, j, k) = 0.;
                }
                else if(k == uCellK_ || cellStatus_(i, j, k + 1) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aT_(i, j, k)*alphaField(i, j, k + 1);
                    aT_(i, j, k) = 0.;
                }

                if(k == 0 && uField.getBottomBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aB_(i, j, k);
                    aB_(i, j, k) = 0.;
                }
                else if(k == 0 || cellStatus_(i, j, k - 1) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aB_(i, j, k)*alphaField(i, j, k - 1);
                    aB_(i, j, k) = 0.;
                }
            }
        }
    }

    // Set-up the solution matrix
    indexMap.generateMap(cellStatus_);
    A.setSize(indexMap.nActive(), indexMap.nActive());
    x.setSize(indexMap.nActive());
    b.setSize(indexMap.nActive());

    A.preallocate(7, 0);

    // Add coefficients to the solution matrix (this process is a little memory intensive, may need to redo this later)
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                A.setValue(indexMap(i, j, k, 0), indexMap(i, j, k, 0), aP_(i, j, k), INSERT_VALUES);

                // I-direction coefficients
                if(i < uCellI_ && cellStatus_(i + 1, j, k) == ACTIVE)
                {
                    A.setValue(indexMap(i, j, k, 0), indexMap(i + 1, j, k, 0), aE_(i, j, k), INSERT_VALUES);
                }

                if(i > 0 && cellStatus_(i - 1, j, k) == ACTIVE)
                {
                    A.setValue(indexMap(i, j, k, 0), indexMap(i - 1, j, k, 0), aW_(i, j, k), INSERT_VALUES);
                }

                // J-direction coefficients
                if(j < uCellJ_ && cellStatus_(i, j + 1, k) == ACTIVE)
                {
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j + 1, k, 0), aN_(i, j, k), INSERT_VALUES);
                }

                if(j > 0 && cellStatus_(i, j - 1, k) == ACTIVE)
                {
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j - 1, k, 0), aS_(i, j, k), INSERT_VALUES);
                }

                // K-direction coefficents
                if(k < uCellK_ && cellStatus_(i, j, k + 1) == ACTIVE)
                {
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j, k + 1, 0), aT_(i, j, k), INSERT_VALUES);
                }

                if(k > 0 && cellStatus_(i, j, k - 1) == ACTIVE)
                {
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j, k - 1, 0), aB_(i, j, k), INSERT_VALUES);
                }

                x.setValue(indexMap(i, j, k, 0), alphaField(i, j, k), INSERT_VALUES);
                b.setValue(indexMap(i, j, k, 0), bP_(i, j, k).x, INSERT_VALUES);
            }
        }
    }

    alphaGmresIters_ = 0;

    // Assemble and solve
    A.assemble();
    x.assemble();
    b.assemble();
    alphaGmresIters_ += A.solve(b, x);

    // Map solution back to the domain
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                alphaField(i, j, k) = x(indexMap(i, j, k, 0));
            }
        }
    }
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
