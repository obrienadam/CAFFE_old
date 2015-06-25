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

#include <algorithm>
#include <iomanip>

#include "Simple.h"
#include "InputStringProcessing.h"

// ************* Constructors and Destructors *************

Simple::Simple()
    :
      uField0_("u0", PRIMITIVE),
      uFieldStar_("uStar", PRIMITIVE),
      a0P_("a0P", AUXILLARY),
      aP_("aP", AUXILLARY),
      aE_("aE", AUXILLARY),
      aW_("aW", AUXILLARY),
      aN_("aN", AUXILLARY),
      aS_("aS", AUXILLARY),
      aT_("aT", AUXILLARY),
      aB_("aB", AUXILLARY),
      dE_("dE", AUXILLARY),
      dW_("dW", AUXILLARY),
      dN_("dN", AUXILLARY),
      dS_("dS", AUXILLARY),
      dT_("dT", AUXILLARY),
      dB_("dB", AUXILLARY),
      cE_("dE", AUXILLARY),
      cW_("dW", AUXILLARY),
      cN_("dN", AUXILLARY),
      cS_("dS", AUXILLARY),
      cT_("dT", AUXILLARY),
      cB_("dB", AUXILLARY),
      bP_("bP", AUXILLARY),
      pCorr_("pCorr", PRIMITIVE),
      gradPCorr_("gradPCorr", PRIMITIVE),
      hField_("hField", PRIMITIVE),
      dField_("dField", PRIMITIVE),
      gradUField_("gradUField", PRIMITIVE),
      gradPField_("gradPField", PRIMITIVE),
      gradVecField_("gradVecField", AUXILLARY),
      gradScalarField_("gradScalarField", AUXILLARY),
      timeAccurate_(false),
      relaxationFactorMomentum_(0.3),
      relaxationFactorPCorr_(0.1),
      gradReconstructionMethod_(DIVERGENCE_THEOREM),
      maxInnerIters_(20),
      momentumGmresIters_(0),
      pCorrGmresIters_(0)
{

}

// ************* Private Methods *************

void Simple::setConstantFields(Input& input)
{
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& rhoField = *rhoFieldPtr_;
    Field<double>& muField = *muFieldPtr_;
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                dE_(i, j, k) = dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                dW_(i, j, k) = dot(mesh.fAreaNormW(i, j, k), mesh.fAreaNormW(i, j, k))/dot(mesh.fAreaNormW(i, j, k), mesh.rCellW(i, j, k));
                dN_(i, j, k) = dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                dS_(i, j, k) = dot(mesh.fAreaNormS(i, j, k), mesh.fAreaNormS(i, j, k))/dot(mesh.fAreaNormS(i, j, k), mesh.rCellS(i, j, k));
                dT_(i, j, k) = dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                dB_(i, j, k) = dot(mesh.fAreaNormB(i, j, k), mesh.fAreaNormB(i, j, k))/dot(mesh.fAreaNormB(i, j, k), mesh.rCellB(i, j, k));

                cE_(i, j, k) = mesh.fAreaNormE(i, j, k) - mesh.rCellE(i, j, k)*dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                cW_(i, j, k) = mesh.fAreaNormW(i, j, k) - mesh.rCellW(i, j, k)*dot(mesh.fAreaNormW(i, j, k), mesh.fAreaNormW(i, j, k))/dot(mesh.fAreaNormW(i, j, k), mesh.rCellW(i, j, k));
                cN_(i, j, k) = mesh.fAreaNormN(i, j, k) - mesh.rCellN(i, j, k)*dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                cS_(i, j, k) = mesh.fAreaNormS(i, j, k) - mesh.rCellS(i, j, k)*dot(mesh.fAreaNormS(i, j, k), mesh.fAreaNormS(i, j, k))/dot(mesh.fAreaNormS(i, j, k), mesh.rCellS(i, j, k));
                cT_(i, j, k) = mesh.fAreaNormT(i, j, k) - mesh.rCellT(i, j, k)*dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                cB_(i, j, k) = mesh.fAreaNormB(i, j, k) - mesh.rCellB(i, j, k)*dot(mesh.fAreaNormB(i, j, k), mesh.fAreaNormB(i, j, k))/dot(mesh.fAreaNormB(i, j, k), mesh.rCellB(i, j, k));

                rhoField(i, j, k) = input.inputDoubles["rho"];
                muField(i, j, k) = input.inputDoubles["mu"];
            }
        }
    }

    rhoField.setAllBoundaries(ZERO_GRADIENT, 0.,
                              ZERO_GRADIENT, 0.,
                              ZERO_GRADIENT, 0.,
                              ZERO_GRADIENT, 0.,
                              ZERO_GRADIENT, 0.,
                              ZERO_GRADIENT, 0.);

    muField.setAllBoundaries(ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.);

    rhoField.setBoundaryFields();
    muField.setBoundaryFields();
    interpolateInteriorFaces(rhoField, NON_WEIGHTED);
    interpolateInteriorFaces(muField, NON_WEIGHTED);
}

void Simple::storeUField(Field<Vector3D> &uField, Field<Vector3D> &uFieldOld)
{
    int i, j, k;

    for(k = 0; k < nFacesK_; ++k)
    {
        for(j = 0; j < nFacesJ_; ++j)
        {
            for(i = 0; i < nFacesI_; ++i)
            {
                if(j < uFaceJ_ && k < uFaceK_)
                    uFieldOld.faceI(i, j, k) = uField.faceI(i, j, k);

                if(i < uFaceI_ && k < uFaceK_)
                    uFieldOld.faceJ(i, j, k) = uField.faceJ(i, j, k);

                if(i < uFaceI_ && j < uFaceJ_)
                    uFieldOld.faceK(i, j, k) = uField.faceK(i, j, k);
            }
        }
    }

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                uFieldOld(i, j, k) = uField(i, j, k);
            }
        }
    }
}

void Simple::computeMomentum(Field<double>& rhoField, Field<double>& muField, Field<double>& massFlowField, Field<Vector3D> *sFieldPtr, double timeStep, Field<Vector3D>& uField, Field<double>& pField)
{
    using namespace std;

    int i, j, k, l;
    HexaFvmMesh& mesh = *meshPtr_;
    SparseMatrix A;
    SparseVector x, b;

    storeUField(uField, uFieldStar_);

    extrapolateInteriorFaces(uField, gradUField_);
    extrapolateInteriorFaces(pField, gradPField_);

    computeCellCenteredGradients(uField, gradUField_, DIVERGENCE_THEOREM);
    computeCellCenteredGradients(pField, gradPField_, DIVERGENCE_THEOREM);

    interpolateInteriorFaces(gradUField_, VOLUME_WEIGHTED);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                // Time coefficient
                if(timeAccurate_)
                    a0P_(i, j, k) = rhoField(i, j, k)*mesh.cellVol(i, j, k)/timeStep;
                else
                    a0P_(i, j, k) = 0.;

                // Face coefficients
                aE_(i, j, k) =  min(massFlowField.faceE(i, j, k), 0.) - muField.faceE(i, j, k)*dE_(i, j, k);
                aW_(i, j, k) =  -max(massFlowField.faceW(i, j, k), 0.) - muField.faceW(i, j, k)*dW_(i, j, k);
                aN_(i, j, k) =  min(massFlowField.faceN(i, j, k), 0.) - muField.faceN(i, j, k)*dN_(i, j, k);
                aS_(i, j, k) =  -max(massFlowField.faceS(i, j, k), 0.) - muField.faceS(i, j, k)*dS_(i, j, k);
                aT_(i, j, k) =  min(massFlowField.faceT(i, j, k), 0.) - muField.faceT(i, j, k)*dT_(i, j, k);
                aB_(i, j, k) =  -max(massFlowField.faceB(i, j, k), 0.) - muField.faceB(i, j, k)*dB_(i, j, k);

                // Central coefficient
                aP_(i, j, k) = (max(massFlowField.faceE(i, j, k), 0.) + muField.faceE(i, j, k)*dE_(i, j, k)
                                - min(massFlowField.faceW(i, j, k), 0.) + muField.faceW(i, j, k)*dW_(i, j, k)
                                + max(massFlowField.faceN(i, j, k), 0.) + muField.faceN(i, j, k)*dN_(i, j, k)
                                - min(massFlowField.faceS(i, j, k), 0.) + muField.faceS(i, j, k)*dS_(i, j, k)
                                + max(massFlowField.faceT(i, j, k), 0.) + muField.faceT(i, j, k)*dT_(i, j, k)
                                - min(massFlowField.faceB(i, j, k), 0.) + muField.faceB(i, j, k)*dB_(i, j, k)
                                + a0P_(i, j, k))/relaxationFactorMomentum_;

                bP_(i, j, k) = a0P_(i, j, k)*uField0_(i, j, k);

                // Compute the cross-diffusion terms (due to mesh non-orthogonality)
                bP_(i, j, k) += muField.faceE(i, j, k)*dot(gradUField_.faceE(i, j, k), cE_(i, j, k));
                bP_(i, j, k) += muField.faceW(i, j, k)*dot(gradUField_.faceW(i, j, k), cW_(i, j, k));
                bP_(i, j, k) += muField.faceN(i, j, k)*dot(gradUField_.faceN(i, j, k), cN_(i, j, k));
                bP_(i, j, k) += muField.faceS(i, j, k)*dot(gradUField_.faceS(i, j, k), cS_(i, j, k));
                bP_(i, j, k) += muField.faceT(i, j, k)*dot(gradUField_.faceT(i, j, k), cT_(i, j, k));
                bP_(i, j, k) += muField.faceB(i, j, k)*dot(gradUField_.faceB(i, j, k), cB_(i, j, k));

                /*
                // Higher-order convection order terms go here (these need to be limited!)
                if(i < uCellI_)
                    bP_(i, j, k) += -min(massFlowField.faceE(i, j, k), 0.)*dot(gradUField_(i + 1, j, k), mesh.rFaceW(i + 1, j, k));

                if(i > 0)
                    bP_(i, j, k) += max(massFlowField.faceW(i, j, k), 0.)*dot(gradUField_(i - 1, j, k), mesh.rFaceE(i - 1, j, k));

                if(j < uCellJ_)
                    bP_(i, j, k) += -min(massFlowField.faceN(i, j, k), 0.)*dot(gradUField_(i, j + 1, k), mesh.rFaceS(i, j + 1, k));

                if(j > 0)
                    bP_(i, j, k) += max(massFlowField.faceS(i, j, k), 0.)*dot(gradUField_(i, j - 1, k), mesh.rFaceN(i, j - 1, k));

                if(k < uCellK_)
                    bP_(i, j, k) += -min(massFlowField.faceT(i, j, k), 0.)*dot(gradUField_(i, j, k + 1), mesh.rFaceB(i, j, k + 1));

                if(k > 0)
                    bP_(i, j, k) += max(massFlowField.faceB(i, j, k), 0.)*dot(gradUField_(i, j, k - 1), mesh.rFaceT(i, j, k - 1));

                bP_(i, j, k) += (-max(massFlowField.faceE(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceE(i, j, k))
                                 + min(massFlowField.faceW(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceW(i, j, k))
                                 - max(massFlowField.faceN(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceN(i, j, k))
                                 + min(massFlowField.faceS(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceS(i, j, k))
                                 - max(massFlowField.faceT(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceT(i, j, k))
                                 + min(massFlowField.faceB(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceB(i, j, k)))/relaxationFactorMomentum_;
                */

                // Pressure term
                bP_(i, j, k) += -mesh.cellVol(i, j, k)*gradPField_(i, j, k);

                // Relaxation source term
                bP_(i, j, k) += (1. - relaxationFactorMomentum_)*aP_(i, j, k)*uFieldStar_(i, j, k);

                // Additional source terms
                if(sFieldPtr != NULL)
                    bP_(i, j, k) += (*sFieldPtr)(i, j, k);

                //- Boundary conditions
                // I-direction bcs
                if(i == uCellI_ && uField.getEastBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aE_(i, j, k);
                    bP_(i, j, k) += -aE_(i, j, k)*dot(gradUField_(i, j, k), mesh.rFaceE(i, j, k));
                    aE_(i, j, k) = 0.;
                }
                else if(i == uCellI_ || cellStatus_(i + 1, j, k) == INTERPOLATION)
                {
                    bP_(i, j, k) += -aE_(i, j, k)*uField(i + 1, j, k);
                    aE_(i, j, k) = 0.;
                }

                if(i == 0 && uField.getWestBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aW_(i, j, k);
                    bP_(i, j, k) += -aW_(i, j, k)*dot(gradUField_(i, j, k), mesh.rFaceW(i, j, k));
                    aW_(i, j, k) = 0.;
                }
                else if(i == 0 || cellStatus_(i - 1, j, k) == INTERPOLATION)
                {
                    bP_(i, j, k) += -aW_(i, j, k)*uField(i - 1, j, k);
                    aW_(i, j, k) = 0.;
                }

                // J-direction bcs
                if(j == uCellJ_ && uField.getNorthBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aN_(i, j, k);
                    bP_(i, j, k) += -aN_(i, j, k)*dot(gradUField_(i, j, k), mesh.rFaceN(i, j, k));
                    aN_(i, j, k) = 0.;
                }
                else if(j == uCellJ_ || cellStatus_(i, j + 1, k) == INTERPOLATION)
                {
                    bP_(i, j, k) += -aN_(i, j, k)*uField(i, j + 1, k);
                    aN_(i, j, k) = 0.;
                }

                if(j == 0 && uField.getSouthBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aS_(i, j, k);
                    bP_(i, j, k) += -aS_(i, j, k)*dot(gradUField_(i, j, k), mesh.rFaceS(i, j, k));
                    aS_(i, j, k) = 0.;
                }
                else if(j == 0 || cellStatus_(i, j - 1, k) == INTERPOLATION)
                {
                    bP_(i, j, k) += -aS_(i, j, k)*uField(i, j - 1, k);
                    aS_(i, j, k) = 0.;
                }

                // K-direction bcs
                if(k == uCellK_ && uField.getTopBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aT_(i, j, k);
                    bP_(i, j, k) += -aT_(i, j, k)*dot(gradUField_(i, j, k), mesh.rFaceT(i, j, k));
                    aT_(i, j, k) = 0.;
                }
                else if(k == uCellK_ || cellStatus_(i, j, k + 1) == INTERPOLATION)
                {
                    bP_(i, j, k) += -aT_(i, j, k)*uField(i, j, k + 1);
                    aT_(i, j, k) = 0.;
                }

                if(k == 0 && uField.getBottomBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aB_(i, j, k);
                    bP_(i, j, k) += -aB_(i, j, k)*dot(gradUField_(i, j, k), mesh.rFaceB(i, j, k));
                    aB_(i, j, k) = 0.;
                }
                else if(k == 0 || cellStatus_(i, j, k - 1) == INTERPOLATION)
                {
                    bP_(i, j, k) += -aB_(i, j, k)*uField(i, j, k - 1);
                    aB_(i, j, k) = 0.;
                }

                // Update D-field
                dField_(i, j, k) = mesh.cellVol(i, j, k)/aP_(i, j, k);
            }
        }
    }

    dField_.setBoundaryFields();

    momentumResidual_ = computeResidual(uFieldStar_);

    // Set-up the solution matrix
    indexMap.generateMap(cellStatus_);
    A.allocate(3*indexMap.nActive(), 3*indexMap.nActive(), 7);
    x.allocate(3*indexMap.nActive());
    b.allocate(3*indexMap.nActive());

    // Add coefficients to the solution matrix (this process is a little memory intensive, may need to redo this later)
    for(l = 0; l < 3; ++l)
    {
        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(cellStatus_(i, j, k) != ACTIVE)
                        continue;

                    A.setValue(indexMap(i, j, k, l), indexMap(i, j, k, l), aP_(i, j, k));

                    // I-direction coefficients
                    if(i < uCellI_ && cellStatus_(i + 1, j, k) == ACTIVE)
                        A.setValue(indexMap(i, j, k, l), indexMap(i + 1, j, k, l), aE_(i, j, k));

                    if(i > 0 && cellStatus_(i - 1, j, k) == ACTIVE)
                        A.setValue(indexMap(i, j, k, l), indexMap(i - 1, j, k, l), aW_(i, j, k));

                    // J-direction coefficients
                    if(j < uCellJ_ && cellStatus_(i, j + 1, k) == ACTIVE)
                        A.setValue(indexMap(i, j, k, l), indexMap(i, j + 1, k, l), aN_(i, j, k));

                    if(j > 0 && cellStatus_(i, j - 1, k) == ACTIVE)
                        A.setValue(indexMap(i, j, k, l), indexMap(i, j - 1, k, l), aS_(i, j, k));

                    // K-direction coefficents
                    if(k < uCellK_ && cellStatus_(i, j, k + 1) == ACTIVE)
                        A.setValue(indexMap(i, j, k, l), indexMap(i, j, k + 1, l), aT_(i, j, k));

                    if(k > 0 && cellStatus_(i, j, k - 1) == ACTIVE)
                        A.setValue(indexMap(i, j, k, l), indexMap(i, j, k - 1, l), aB_(i, j, k));

                    b.setValue(indexMap(i, j, k, l), bP_(i, j, k)(l));
                }
            }
        }
    }

    momentumGmresIters_ += A.solve(b, x);

    // Map solution back to the domain
    for(l = 0; l < 3; ++l)
    {
        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(cellStatus_(i, j, k) != ACTIVE)
                        continue;

                    uField(i, j, k)(l) = x(indexMap(i, j, k, l));
                }
            }
        }
    }

    // Assemble the pseudo-velocity vector, to be used for the Rhie-Chow interpolation
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                hField_(i, j, k) = (- aE_(i, j, k)*uField(i + 1, j, k)
                                    - aW_(i, j, k)*uField(i - 1, j, k)
                                    - aN_(i, j, k)*uField(i, j + 1, k)
                                    - aS_(i, j, k)*uField(i, j - 1, k)
                                    - aT_(i, j, k)*uField(i, j, k + 1)
                                    - aB_(i, j, k)*uField(i, j, k - 1)
                                    + bP_(i, j, k)
                                    + mesh.cellVol(i, j, k)*gradPField_(i, j, k))/aP_(i, j, k);
            }
        }
    }

    hField_.setBoundaryFields();
}

void Simple::rhieChowInterpolateInteriorFaces(Field<Vector3D> &uField, Field<double>& pField)
{
    int i, j, k;
    Vector3D sf, ds;

    extrapolateInteriorFaces(hField_, gradVecField_);
    extrapolateInteriorFaces(dField_, gradScalarField_);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                // Check if west, south and bottom bcs require the face mass-flux to be computed
                if(i == 0 && uField.getWestBoundaryPatch() == ZERO_GRADIENT)
                {
                    getFaceStencil(i, j, k, WEST, sf, ds);
                    uField.faceW(i, j, k) = hField_.faceW(i, j, k) - dField_.faceW(i, j, k)*(pField(i - 1, j, k) - pField(i, j, k))*sf/dot(sf, ds);
                }

                if(j == 0 && uField.getSouthBoundaryPatch() == ZERO_GRADIENT)
                {
                    getFaceStencil(i, j, k, SOUTH, sf, ds);
                    uField.faceS(i, j, k) = hField_.faceS(i, j, k) - dField_.faceS(i, j, k)*(pField(i, j - 1, k) - pField(i, j, k))*sf/dot(sf, ds);
                }

                if(k == 0 && uField.getBottomBoundaryPatch() == ZERO_GRADIENT)
                {
                    getFaceStencil(i, j, k, BOTTOM, sf, ds);
                    uField.faceB(i, j, k) = hField_.faceB(i, j, k) - dField_.faceB(i, j, k)*(pField(i, j, k - 1) - pField(i, j, k))*sf/dot(sf, ds);
                }

                // Check if the east, south and top bcs require the face mass-flux to be computed
                if(i < uCellI_ || uField.getEastBoundaryPatch() == ZERO_GRADIENT)
                {
                    getFaceStencil(i, j, k, EAST, sf, ds);
                    uField.faceE(i, j, k) = hField_.faceE(i, j, k) - dField_.faceE(i, j, k)*(pField(i + 1, j, k) - pField(i, j, k))*sf/dot(sf, ds);
                }

                if(j < uCellJ_ || uField.getNorthBoundaryPatch() == ZERO_GRADIENT)
                {
                    getFaceStencil(i, j, k, NORTH, sf, ds);
                    uField.faceN(i, j, k) = hField_.faceN(i, j, k) - dField_.faceN(i, j, k)*(pField(i, j + 1, k) - pField(i, j, k))*sf/dot(sf, ds);
                }

                if(k < uCellK_ || uField.getTopBoundaryPatch() == ZERO_GRADIENT)
                {
                    getFaceStencil(i, j, k, TOP, sf, ds);
                    uField.faceT(i, j, k) = hField_.faceT(i, j, k) - dField_.faceT(i, j, k)*(pField(i, j, k + 1) - pField(i, j, k))*sf/dot(sf, ds);
                }
            }
        }
    }
}

void Simple::computeMassFlowFaces(Field<double>& rhoField, Field<Vector3D> &uField, Field<double>& massFlowField)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;

    for(k = 0; k < nFacesK_; ++k)
    {
        for(j = 0; j < nFacesJ_; ++j)
        {
            for(i = 0; i < nFacesI_; ++i)
            {
                if(j < uFaceJ_ && k < uFaceK_)
                    massFlowField.faceI(i, j, k) = rhoField.faceI(i, j, k)*dot(uField.faceI(i, j, k), mesh.fAreaNormI(i, j, k));

                if(i < uFaceI_ && k < uFaceK_)
                    massFlowField.faceJ(i, j, k) = rhoField.faceJ(i, j, k)*dot(uField.faceJ(i, j, k), mesh.fAreaNormJ(i, j, k));

                if(i < uFaceI_ && j < uFaceJ_)
                    massFlowField.faceK(i, j, k) = rhoField.faceK(i, j, k)*dot(uField.faceK(i, j, k), mesh.fAreaNormK(i, j, k));
            }
        }
    }
}

void Simple::computePCorr(Field<double>& rhoField, Field<double>& massFlowField, Field<Vector3D>& uField, Field<double>& pField)
{
    int i, j, k;
    Vector3D sf, ds, gradPCorrBar;
    SparseMatrix A;
    SparseVector x, b;

    rhieChowInterpolateInteriorFaces(uField, pField);
    computeMassFlowFaces(rhoField, uField, massFlowField);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                aE_(i, j, k) = rhoField.faceE(i, j, k)*dField_.faceE(i, j, k)*dE_(i, j, k);
                aW_(i, j, k) = rhoField.faceW(i, j, k)*dField_.faceW(i, j, k)*dW_(i, j, k);
                aN_(i, j, k) = rhoField.faceN(i, j, k)*dField_.faceN(i, j, k)*dN_(i, j, k);
                aS_(i, j, k) = rhoField.faceS(i, j, k)*dField_.faceS(i, j, k)*dS_(i, j, k);
                aT_(i, j, k) = rhoField.faceT(i, j, k)*dField_.faceT(i, j, k)*dT_(i, j, k);
                aB_(i, j, k) = rhoField.faceB(i, j, k)*dField_.faceB(i, j, k)*dB_(i, j, k);

                aP_(i, j, k) = -(aE_(i, j, k) + aW_(i, j, k) + aN_(i, j, k) + aS_(i, j, k) + aT_(i, j, k) + aB_(i, j, k));

                massFlowField(i, j, k) = massFlowField.faceE(i, j, k) - massFlowField.faceW(i, j, k)
                        + massFlowField.faceN(i, j, k) - massFlowField.faceS(i, j, k)
                        + massFlowField.faceT(i, j, k) - massFlowField.faceB(i, j, k);

                bP_(i, j, k).x = massFlowField(i, j, k);

                //- Boundary conditions
                // I-direction bcs
                if(i == uCellI_ && pCorr_.getEastBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aE_(i, j, k);
                    aE_(i, j, k) = 0.;
                }
                else if(i == uCellI_ || cellStatus_(i + 1, j, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aE_(i, j, k)*pCorr_(i + 1, j, k);
                    aE_(i, j, k) = 0.;
                }

                if(i == 0 && pCorr_.getWestBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aW_(i, j, k);
                    aW_(i, j, k) = 0.;
                }
                else if(i == 0 || cellStatus_(i - 1, j, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aW_(i, j, k)*pCorr_(i - 1, j, k);
                    aW_(i, j, k) = 0.;
                }

                // J-direction bcs
                if(j == uCellJ_ && pCorr_.getNorthBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aN_(i, j, k);
                    aN_(i, j, k) = 0.;
                }
                else if(j == uCellJ_ || cellStatus_(i, j + 1, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aN_(i, j, k)*pCorr_(i, j + 1, k);
                    aN_(i, j, k) = 0.;
                }

                if(j == 0 && pCorr_.getSouthBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aS_(i, j, k);
                    aS_(i, j, k) = 0.;
                }
                else if(j == 0 || cellStatus_(i, j - 1, k) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aS_(i, j, k)*pCorr_(i, j - 1, k);
                    aS_(i, j, k) = 0.;
                }

                // K-direction bcs
                if(k == uCellK_ && pCorr_.getTopBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aT_(i, j, k);
                    aT_(i, j, k) = 0.;
                }
                else if(k == uCellK_ || cellStatus_(i, j, k + 1) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aT_(i, j, k)*pCorr_(i, j, k + 1);
                    aT_(i, j, k) = 0.;
                }

                if(k == 0 && pCorr_.getBottomBoundaryPatch() == ZERO_GRADIENT)
                {
                    aP_(i, j, k) += aB_(i, j, k);
                    aB_(i, j, k) = 0.;
                }
                else if(k == 0 || cellStatus_(i, j, k - 1) == INTERPOLATION)
                {
                    bP_(i, j, k).x += -aB_(i, j, k)*pCorr_(i, j, k - 1);
                    aB_(i, j, k) = 0.;
                }
            }
        }
    }

    // Set-up the solution matrix
    indexMap.generateMap(cellStatus_);
    A.allocate(indexMap.nActive(), indexMap.nActive(), 7);
    x.allocate(indexMap.nActive());
    b.allocate(indexMap.nActive());

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                A.setValue(indexMap(i, j, k, 0), indexMap(i, j, k, 0), aP_(i, j, k));

                // I-direction coefficients
                if(i < uCellI_ && cellStatus_(i + 1, j, k) == ACTIVE)
                    A.setValue(indexMap(i, j, k, 0), indexMap(i + 1, j, k, 0), aE_(i, j, k));

                if(i > 0 && cellStatus_(i - 1, j, k) == ACTIVE)
                    A.setValue(indexMap(i, j, k, 0), indexMap(i - 1, j, k, 0), aW_(i, j, k));

                // J-direction coefficients
                if(j < uCellJ_ && cellStatus_(i, j + 1, k) == ACTIVE)
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j + 1, k, 0), aN_(i, j, k));

                if(j > 0 && cellStatus_(i, j - 1, k) == ACTIVE)
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j - 1, k, 0), aS_(i, j, k));

                // K-direction coefficients
                if(k < uCellK_ && cellStatus_(i, j, k + 1) == ACTIVE)
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j, k + 1, 0), aT_(i, j, k));

                if(k > 0 && cellStatus_(i, j, k - 1) == ACTIVE)
                    A.setValue(indexMap(i, j, k, 0), indexMap(i, j, k - 1, 0), aB_(i, j, k));

                b.setValue(indexMap(i, j, k, 0), bP_(i, j, k).x);
            }
        }
    }

    pCorrGmresIters_ += A.solve(b, x);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                pCorr_(i, j, k) = x(indexMap(i, j, k, 0));
            }
        }
    }
}

void Simple::correctContinuity(Field<double>& rhoField, Field<double>& massFlowField, Field<Vector3D> &uField, Field<double> &pField)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& massFlow = mesh.findScalarField("massFlow");

    extrapolateInteriorFaces(pCorr_, gradPCorr_);
    computeCellCenteredGradients(pCorr_, gradPCorr_, DIVERGENCE_THEOREM);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                //- Correct mass flow
                // Correct west, south and bottom boundary mass flows if they are outlets
                if(i == 0 && uField.getWestBoundaryPatch() == ZERO_GRADIENT)
                    massFlowField.faceW(i, j, k) += -rhoField.faceW(i, j, k)*dField_.faceW(i, j, k)*dW_(i, j, k)*(pCorr_(i - 1, j, k) - pCorr_(i, j, k));

                if(j == 0 && uField.getSouthBoundaryPatch() == ZERO_GRADIENT)
                    massFlowField.faceS(i, j, k) += -rhoField.faceS(i, j, k)*dField_.faceS(i, j, k)*dS_(i, j, k)*(pCorr_(i, j - 1, k) - pCorr_(i, j, k));

                if(k == 0 && uField.getBottomBoundaryPatch() == ZERO_GRADIENT)
                    massFlowField.faceB(i, j, k) += -rhoField.faceB(i, j, k)*dField_.faceB(i, j, k)*dB_(i, j, k)*(pCorr_(i, j, k - 1) - pCorr_(i, j, k));

                // Correct all interior faces and east, north and top boundary mass flows if they are outlets
                if(i < uCellI_ || uField.getEastBoundaryPatch() == ZERO_GRADIENT)
                    massFlowField.faceE(i, j, k) += -rhoField.faceE(i, j, k)*dField_.faceE(i, j, k)*dE_(i, j, k)*(pCorr_(i + 1, j, k) - pCorr_(i, j, k));

                if(j < uCellJ_ || uField.getNorthBoundaryPatch() == ZERO_GRADIENT)
                    massFlowField.faceN(i, j, k) += -rhoField.faceN(i, j, k)*dField_.faceN(i, j, k)*dN_(i, j, k)*(pCorr_(i, j + 1, k) - pCorr_(i, j, k));

                if(k < uCellK_ || uField.getTopBoundaryPatch() == ZERO_GRADIENT)
                    massFlowField.faceT(i, j, k) += -rhoField.faceT(i, j, k)*dField_.faceT(i, j, k)*dT_(i, j, k)*(pCorr_(i, j, k + 1) - pCorr_(i, j, k));

                massFlow(i, j, k) = massFlowField.faceE(i, j, k) - massFlowField.faceW(i, j, k)
                        + massFlowField.faceN(i, j, k) - massFlowField.faceS(i, j, k)
                        + massFlowField.faceT(i, j, k) - massFlowField.faceB(i, j, k);

                massFlow(i, j, k) = fabs(massFlow(i, j, k));

                // Correct the pressure field
                pField(i, j, k) += relaxationFactorPCorr_*pCorr_(i, j, k);

                // Correct the velocity field
                uField(i, j, k) += -dField_(i, j, k)*gradPCorr_(i, j, k);
            }
        }
    }

    pField.setBoundaryFields();
    uField.setBoundaryFields();
}

Vector3D Simple::computeResidual(Field<Vector3D> &uField)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;
    Vector3D sumMomentumResidualSqr = Vector3D(0., 0., 0.);
    double sumVol = 0.;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                sumVol += mesh.cellVol(i, j, k);
                sumMomentumResidualSqr += mesh.cellVol(i, j, k)*sqr(aP_(i, j, k)*uField(i, j, k)
                                                                    + aE_(i, j, k)*uField(i + 1, j, k)
                                                                    + aW_(i, j, k)*uField(i - 1, j, k)
                                                                    + aN_(i, j, k)*uField(i, j + 1, k)
                                                                    + aS_(i, j, k)*uField(i, j - 1, k)
                                                                    + aT_(i, j, k)*uField(i, j, k + 1)
                                                                    + aB_(i, j, k)*uField(i, j, k - 1)
                                                                    - bP_(i, j, k));
            }
        }
    }

    if(isnan(sumMomentumResidualSqr.x) || isnan(sumMomentumResidualSqr.y) || isnan(sumMomentumResidualSqr.z))
        Output::raiseException("Simple", "computeResidual", "a NaN value was detected.");

    return sqrt(sumMomentumResidualSqr/sumVol);
}

// ************* Public Methods *************

void Simple::initialize(Input &input, HexaFvmMesh &mesh)
{
    FvScheme::initialize(input, mesh, "NA");

    //- Set respective index parameters
    uxStartI_ = 0;
    uxEndI_ = nCells_ - 1;

    uyStartI_ = uxEndI_ + 1;
    uyEndI_ = uyStartI_ + nCells_ - 1;

    uzStartI_ = uyEndI_ + 1;
    uzEndI_ = uzStartI_ + nCells_ - 1;

    pStartI_ = uzEndI_ + 1;
    pEndI_ = pStartI_ + nCells_ - 1;

    //- Find the relevant fields
    uFieldPtr_ = &mesh.findVectorField("u");
    pFieldPtr_ = &mesh.findScalarField("p");
    rhoFieldPtr_ = &mesh.findScalarField("rho");
    muFieldPtr_ = &mesh.findScalarField("mu");
    massFlowFieldPtr_ = &mesh.findScalarField("massFlow");

    //- Read all relevant input parameters
    if(input.inputStrings["timeAccurate"] == "ON")
        timeAccurate_ = true;
    else if (input.inputStrings["timeAccurate"] == "OFF")
        timeAccurate_ = false;

    relaxationFactorMomentum_ = input.inputDoubles["relaxationFactorMomentum"];
    relaxationFactorPCorr_ = input.inputDoubles["relaxationFactorPCorr"];

    maxInnerIters_ = input.inputInts["maxInnerIters"];

    //- Allocate memory
    uField0_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    uFieldStar_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    a0P_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aP_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aE_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aW_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aN_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aS_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aT_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aB_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    dE_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dW_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dN_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dS_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dT_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dB_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    cE_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cW_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cN_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cS_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cT_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cB_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    bP_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    hField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dField_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    pCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    gradPCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradUField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradPField_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    gradVecField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradScalarField_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    //- Set the boundary conditions
    setBoundaryConditions(input);

    //- Set constant fields
    setConstantFields(input);

    Output::print("Simple", "Time accurate : " + input.inputStrings["timeAccurate"]);
    Output::print("Simple", "initialization complete.");
}

void Simple::setBoundaryConditions(Input &input)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    BoundaryPatch uBoundaryTypes[6], pBoundaryTypes[6];
    Vector3D uBoundaryRefValues[6];
    std::string key1, key2, flowBoundaryType, buffer;

    for(int i = 0; i < 6; ++i)
    {
        switch(i)
        {
        case 0:
            key1 = "boundaryTypeEast";
            key2 = "boundaryRefValueEast";
            break;
        case 1:
            key1 = "boundaryTypeWest";
            key2 = "boundaryRefValueWest";
            break;
        case 2:
            key1 = "boundaryTypeNorth";
            key2 = "boundaryRefValueNorth";
            break;
        case 3:
            key1 = "boundaryTypeSouth";
            key2 = "boundaryRefValueSouth";
            break;
        case 4:
            key1 = "boundaryTypeTop";
            key2 = "boundaryRefValueTop";
            break;
        case 5:
            key1 = "boundaryTypeBottom";
            key2 = "boundaryRefValueBottom";
            break;
        };

        flowBoundaryType = input.inputStrings[key1];

        buffer = input.inputStrings[key2];
        InputStringProcessing::processBuffer(buffer, true);

        uBoundaryRefValues[i] = stov(buffer);

        if(flowBoundaryType == "inlet")
        {
            uBoundaryTypes[i] = FIXED;
            pBoundaryTypes[i] = ZERO_GRADIENT;
        }
        else if(flowBoundaryType == "outlet")
        {
            uBoundaryTypes[i] = ZERO_GRADIENT;
            pBoundaryTypes[i] = FIXED;
        }
        else if (flowBoundaryType == "wall")
        {
            uBoundaryTypes[i] = FIXED;
            pBoundaryTypes[i] = ZERO_GRADIENT;
        }
        else
        {
            Output::raiseException("Simple", "setBoundaryConditions", "boundary type \"" + flowBoundaryType + "\" is not a valid flow boundary.");
        }
    }

    uField.setAllBoundaries(uBoundaryTypes[0], uBoundaryRefValues[0],
            uBoundaryTypes[1], uBoundaryRefValues[1],
            uBoundaryTypes[2], uBoundaryRefValues[2],
            uBoundaryTypes[3], uBoundaryRefValues[3],
            uBoundaryTypes[4], uBoundaryRefValues[4],
            uBoundaryTypes[5], uBoundaryRefValues[5]);

    pField.setAllBoundaries(pBoundaryTypes[0], 0.,
            pBoundaryTypes[1], 0.,
            pBoundaryTypes[2], 0.,
            pBoundaryTypes[3], 0.,
            pBoundaryTypes[4], 0.,
            pBoundaryTypes[5], 0.);

    pCorr_.setAllBoundaries(pBoundaryTypes[0], 0.,
            pBoundaryTypes[1], 0.,
            pBoundaryTypes[2], 0.,
            pBoundaryTypes[3], 0.,
            pBoundaryTypes[4], 0.,
            pBoundaryTypes[5], 0.);

    dField_.setAllBoundaries(ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.);

    hField_.setAllBoundaries(ZERO_GRADIENT, Vector3D(0., 0., 0.),
                             ZERO_GRADIENT, Vector3D(0., 0., 0.),
                             ZERO_GRADIENT, Vector3D(0., 0., 0.),
                             ZERO_GRADIENT, Vector3D(0., 0., 0.),
                             ZERO_GRADIENT, Vector3D(0., 0., 0.),
                             ZERO_GRADIENT, Vector3D(0., 0., 0.));
}

int Simple::nConservedVariables()
{
    return 4*meshPtr_->size();
}

void Simple::discretize(double timeStep, std::vector<double> &timeDerivatives)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    Field<double>& rhoField = *rhoFieldPtr_;
    Field<double>& muField = *muFieldPtr_;
    Field<double>& massFlowField = *massFlowFieldPtr_;
    int i;

    storeUField(uField, uField0_);

    for(i = 0; i < maxInnerIters_; ++i)
    {
        computeMomentum(rhoField, muField, massFlowField, NULL, timeStep, uField, pField);
        computePCorr(rhoField, massFlowField, uField, pField);
        correctContinuity(rhoField, massFlowField, uField, pField);
        std::cout << "\rDTS iteration completion  |      " << (i + 1) << "/" << maxInnerIters_ << std::fixed << std::setprecision(2) << " (" << 100.*(i + 1)/maxInnerIters_ << "%)";
    }
}

void Simple::copySolution(std::vector<double> &original)
{

}

void Simple::updateSolution(std::vector<double> &update, int method)
{

}

void Simple::displayUpdateMessage()
{
    int i, j, k;
    double maxContinuityError = 0.;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                maxContinuityError = 0.;
            }
        }
    }

    Output::print("Simple", "Momentum prediction total BiCGStab iterations : " + std::to_string(momentumGmresIters_));
    Output::print("Simple", "Pressure correction total BiCGStab iterations : " + std::to_string(pCorrGmresIters_) + "\n");
    Output::print("Simple", "U-Momentum residual           : " + std::to_string(momentumResidual_.x));
    Output::print("Simple", "V-Momentum residual           : " + std::to_string(momentumResidual_.y));
    Output::print("Simple", "W-Momentum residual           : " + std::to_string(momentumResidual_.z));
    Output::print("Simple", "Maximum cell continuity error : " + std::to_string(maxContinuityError) + "\n");

    momentumGmresIters_ = 0;
    pCorrGmresIters_ = 0;
}
