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
      momentumBiCGStabIters_(0),
      pCorrBiCGStabIters_(0)
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

void Simple::createMatrices()
{
    indexMap.generateMap(cellStatus_);

    AMomentum.allocate(3*indexMap.nActive(), 3*indexMap.nActive(), 7);
    xMomentum.allocate(3*indexMap.nActive());
    bMomentum.allocate(3*indexMap.nActive());
    rMomentum.allocate(3*indexMap.nActive());

    APCorr.allocate(indexMap.nActive(), indexMap.nActive(), 7);
    xPCorr.allocate(indexMap.nActive());
    bPCorr.allocate(indexMap.nActive());
}

void Simple::destroyMatrices()
{
    AMomentum.deallocate();
    xMomentum.deallocate();
    bMomentum.deallocate();
    rMomentum.deallocate();

    APCorr.deallocate();
    xPCorr.deallocate();
    bPCorr.deallocate();
}

void Simple::zeroMatrices()
{
    AMomentum.zeroEntries();
    xMomentum.zeroEntries();
    bMomentum.zeroEntries();
    rMomentum.zeroEntries();

    APCorr.zeroEntries();
    xPCorr.zeroEntries();
    bPCorr.zeroEntries();
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
    double a0P, aP, aE, aW, aN, aS, aT, aB;
    double a[7];
    int rowNo, cols[7];
    Vector3D bP(0., 0., 0.);

    storeUField(uField, uFieldStar_);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                // Time coefficient
                if(timeAccurate_)
                    a0P = rhoField(i, j, k)*mesh.cellVol(i, j, k)/timeStep;
                else
                    a0P = 0.;

                // Face coefficients
                aE =  min(massFlowField.faceE(i, j, k), 0.) - muField.faceE(i, j, k)*dE_(i, j, k);
                aW = -max(massFlowField.faceW(i, j, k), 0.) - muField.faceW(i, j, k)*dW_(i, j, k);
                aN =  min(massFlowField.faceN(i, j, k), 0.) - muField.faceN(i, j, k)*dN_(i, j, k);
                aS = -max(massFlowField.faceS(i, j, k), 0.) - muField.faceS(i, j, k)*dS_(i, j, k);
                aT =  min(massFlowField.faceT(i, j, k), 0.) - muField.faceT(i, j, k)*dT_(i, j, k);
                aB = -max(massFlowField.faceB(i, j, k), 0.) - muField.faceB(i, j, k)*dB_(i, j, k);

                // Central coefficient
                aP = (max(massFlowField.faceE(i, j, k), 0.) + muField.faceE(i, j, k)*dE_(i, j, k)
                                - min(massFlowField.faceW(i, j, k), 0.) + muField.faceW(i, j, k)*dW_(i, j, k)
                                + max(massFlowField.faceN(i, j, k), 0.) + muField.faceN(i, j, k)*dN_(i, j, k)
                                - min(massFlowField.faceS(i, j, k), 0.) + muField.faceS(i, j, k)*dS_(i, j, k)
                                + max(massFlowField.faceT(i, j, k), 0.) + muField.faceT(i, j, k)*dT_(i, j, k)
                                - min(massFlowField.faceB(i, j, k), 0.) + muField.faceB(i, j, k)*dB_(i, j, k)
                                + a0P)/relaxationFactorMomentum_;

                bP = a0P*uField0_(i, j, k);

                // Compute the cross-diffusion terms (due to mesh non-orthogonality)
                bP += muField.faceE(i, j, k)*dot(gradUField_.faceE(i, j, k), cE_(i, j, k));
                bP += muField.faceW(i, j, k)*dot(gradUField_.faceW(i, j, k), cW_(i, j, k));
                bP += muField.faceN(i, j, k)*dot(gradUField_.faceN(i, j, k), cN_(i, j, k));
                bP += muField.faceS(i, j, k)*dot(gradUField_.faceS(i, j, k), cS_(i, j, k));
                bP += muField.faceT(i, j, k)*dot(gradUField_.faceT(i, j, k), cT_(i, j, k));
                bP += muField.faceB(i, j, k)*dot(gradUField_.faceB(i, j, k), cB_(i, j, k));

                /*
                // Higher-order convection order terms go here (these need to be limited!)
                if(i < uCellI_)
                    bP += -min(massFlowField.faceE(i, j, k), 0.)*dot(gradUField_(i + 1, j, k), mesh.rFaceW(i + 1, j, k));

                if(i > 0)
                    bP += max(massFlowField.faceW(i, j, k), 0.)*dot(gradUField_(i - 1, j, k), mesh.rFaceE(i - 1, j, k));

                if(j < uCellJ_)
                    bP += -min(massFlowField.faceN(i, j, k), 0.)*dot(gradUField_(i, j + 1, k), mesh.rFaceS(i, j + 1, k));

                if(j > 0)
                    bP += max(massFlowField.faceS(i, j, k), 0.)*dot(gradUField_(i, j - 1, k), mesh.rFaceN(i, j - 1, k));

                if(k < uCellK_)
                    bP += -min(massFlowField.faceT(i, j, k), 0.)*dot(gradUField_(i, j, k + 1), mesh.rFaceB(i, j, k + 1));

                if(k > 0)
                    bP += max(massFlowField.faceB(i, j, k), 0.)*dot(gradUField_(i, j, k - 1), mesh.rFaceT(i, j, k - 1));

                bP += (-max(massFlowField.faceE(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceE(i, j, k))
                                 + min(massFlowField.faceW(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceW(i, j, k))
                                 - max(massFlowField.faceN(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceN(i, j, k))
                                 + min(massFlowField.faceS(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceS(i, j, k))
                                 - max(massFlowField.faceT(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceT(i, j, k))
                                 + min(massFlowField.faceB(i, j, k), 0.)*dot(gradUField_(i, j, k), mesh.rFaceB(i, j, k)))/relaxationFactorMomentum_;
                */

                // Pressure term
                bP += -mesh.cellVol(i, j, k)*gradPField_(i, j, k);

                // Relaxation source term
                bP += (1. - relaxationFactorMomentum_)*aP*uFieldStar_(i, j, k);

                // Additional source terms
                if(sFieldPtr != NULL)
                    bP += (*sFieldPtr)(i, j, k);

                setFieldBcCoeffs(uField, i, j, k, aP, aE, aW, aN, aS, aT, aB, bP);

                a[0] = aP;
                a[1] = aE;
                a[2] = aW;
                a[3] = aN;
                a[4] = aS;
                a[5] = aT;
                a[6] = aB;

                for(l = 0; l < 3; ++l)
                {
                    rowNo = indexMap(i, j, k, l);
                    cols[0] = indexMap(i, j, k, l);
                    cols[1] = indexMap(i + 1, j, k, l);
                    cols[2] = indexMap(i - 1, j, k, l);
                    cols[3] = indexMap(i, j + 1, k, l);
                    cols[4] = indexMap(i, j - 1, k, l);
                    cols[5] = indexMap(i, j, k + 1, l);
                    cols[6] = indexMap(i, j, k - 1, l);
                    AMomentum.setRow(rowNo, 7, cols, a);
                    bMomentum.setValue(rowNo, bP(l));
                }

                // Update D-field
                dField_(i, j, k) = mesh.cellVol(i, j, k)/aP;
            }
        }
    }

    momentumBiCGStabIters_ += AMomentum.solve(bMomentum, xMomentum);

    // Map solution back to the domain
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                uField(i, j, k).x = xMomentum(indexMap(i, j, k, 0));
                uField(i, j, k).y = xMomentum(indexMap(i, j, k, 1));
                uField(i, j, k).z = xMomentum(indexMap(i, j, k, 2));
            }
        }
    }

    // Compute the pseudo-velocity vector
    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                hField_(i, j, k) = uField(i, j, k) + dField_(i, j, k)*gradPField_(i, j, k);
            }
        }
    }

    dField_.setBoundaryFields();
    extrapolateInteriorFaces(dField_, gradScalarField_);

    hField_.setBoundaryFields();
    extrapolateInteriorFaces(hField_, gradVecField_);

    uField.setBoundaryFields();
    rhieChowInterpolateInteriorFaces(uField, pField);
    computeMassFlowFaces(rhoField, uField, massFlowField);
}

void Simple::rhieChowInterpolateInteriorFaces(Field<Vector3D> &uField, Field<double>& pField)
{
    int i, j, k;
    Vector3D sf, ds;

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
    double aP, aE, aW, aN, aS, aT, aB, bP = 0., a[7];
    int rowNo, cols[7];

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                aE = rhoField.faceE(i, j, k)*dField_.faceE(i, j, k)*dE_(i, j, k);
                aW = rhoField.faceW(i, j, k)*dField_.faceW(i, j, k)*dW_(i, j, k);
                aN = rhoField.faceN(i, j, k)*dField_.faceN(i, j, k)*dN_(i, j, k);
                aS = rhoField.faceS(i, j, k)*dField_.faceS(i, j, k)*dS_(i, j, k);
                aT = rhoField.faceT(i, j, k)*dField_.faceT(i, j, k)*dT_(i, j, k);
                aB = rhoField.faceB(i, j, k)*dField_.faceB(i, j, k)*dB_(i, j, k);

                aP = -(aE + aW + aN + aS + aT + aB);

                bP = massFlowField.faceE(i, j, k) - massFlowField.faceW(i, j, k)
                        + massFlowField.faceN(i, j, k) - massFlowField.faceS(i, j, k)
                        + massFlowField.faceT(i, j, k) - massFlowField.faceB(i, j, k);

                setFieldBcCoeffs(pCorr_, i, j, k, aP, aE, aW, aN, aS, aT, aB, bP);

                a[0] = aP;
                a[1] = aE;
                a[2] = aW;
                a[3] = aN;
                a[4] = aS;
                a[5] = aT;
                a[6] = aB;

                rowNo = indexMap(i, j, k, 0);
                cols[0] = indexMap(i, j, k, 0);
                cols[1] = indexMap(i + 1, j, k, 0);
                cols[2] = indexMap(i - 1, j, k, 0);
                cols[3] = indexMap(i, j + 1, k, 0);
                cols[4] = indexMap(i, j - 1, k, 0);
                cols[5] = indexMap(i, j, k + 1, 0);
                cols[6] = indexMap(i, j, k - 1, 0);
                APCorr.setRow(rowNo, 7, cols, a);
                bPCorr.setValue(rowNo, bP);
            }
        }
    }

    pCorrBiCGStabIters_ += APCorr.solve(bPCorr, xPCorr);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                pCorr_(i, j, k) = xPCorr(indexMap(i, j, k, 0));
            }
        }
    }

    pCorr_.setBoundaryFields();
    extrapolateInteriorFaces(pCorr_, gradPCorr_);
    computeCellCenteredGradients(pCorr_, gradPCorr_, DIVERGENCE_THEOREM);
}

void Simple::correctContinuity(Field<double>& rhoField, Field<double>& massFlowField, Field<Vector3D> &uField, Field<double> &pField)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& massFlowError = mesh.findScalarField("massFlowError");

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

                massFlowError(i, j, k) = fabs(massFlowField.faceE(i, j, k) - massFlowField.faceW(i, j, k)
                                              + massFlowField.faceN(i, j, k) - massFlowField.faceS(i, j, k)
                                              + massFlowField.faceT(i, j, k) - massFlowField.faceB(i, j, k));

                // Correct the pressure field
                pField(i, j, k) += relaxationFactorPCorr_*pCorr_(i, j, k);

                // Correct the velocity field
                uField(i, j, k) += -dField_(i, j, k)*gradPCorr_(i, j, k);
            }
        }
    }

    pField.setBoundaryFields();
    extrapolateInteriorFaces(pField, gradPField_);
    computeCellCenteredGradients(pField, gradPField_, DIVERGENCE_THEOREM);

    uField.setBoundaryFields();
    extrapolateInteriorFaces(uField, gradUField_);
    computeCellCenteredGradients(uField, gradUField_, DIVERGENCE_THEOREM);
}

void Simple::computeResidual(Field<Vector3D> &uField)
{
    int i, j, k, l;

    // Momentum residuals

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                for(l = 0; l < 3; ++l)
                {
                    xMomentum.setValue(indexMap(i, j, k, l), uField(i, j, k)(l));
                }
            }
        }
    }

    scale(-1., bMomentum);
    multiplyAdd(AMomentum, xMomentum, bMomentum, rMomentum);
    momentumL1Norm_ = rMomentum.l1Norm();
    momentumL2Norm_ = rMomentum.l2Norm();
    momentumInfNorm_ = rMomentum.infNorm();
}

// ************* Public Methods *************

void Simple::initialize(Input &input, HexaFvmMesh &mesh)
{
    FvScheme::initialize(input, mesh, "NA");

    //- Find the relevant fields
    uFieldPtr_ = &mesh.findVectorField("u");
    pFieldPtr_ = &mesh.findScalarField("p");
    rhoFieldPtr_ = &mesh.findScalarField("rho");
    muFieldPtr_ = &mesh.findScalarField("mu");
    massFlowFieldPtr_ = &mesh.findScalarField("massFlowError");

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

    //- Create matrices
    createMatrices();

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
        else if(flowBoundaryType == "empty")
        {
            uBoundaryTypes[i] = EMPTY;
            pBoundaryTypes[i] = EMPTY;
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
        zeroMatrices();
        computeMomentum(rhoField, muField, massFlowField, NULL, timeStep, uField, pField);
        computePCorr(rhoField, massFlowField, uField, pField);
        correctContinuity(rhoField, massFlowField, uField, pField);
        std::cout << "\rDTS iteration completion  |      " << (i + 1) << "/" << maxInnerIters_ << std::fixed << std::setprecision(2) << " (" << 100.*(i + 1)/maxInnerIters_ << "%)";
    }

    computeResidual(uField);
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

    Output::print("Simple", "Momentum prediction total BiCGStab iterations : " + std::to_string(momentumBiCGStabIters_));
    Output::print("Simple", "Pressure correction total BiCGStab iterations : " + std::to_string(pCorrBiCGStabIters_) + "\n");
    Output::print("Simple", "Momentum L1 residual norm     : " + std::to_string(momentumL1Norm_));
    Output::print("Simple", "Momentum L2 residual norm     : " + std::to_string(momentumL2Norm_));
    Output::print("Simple", "Momentum inf residual norm    : " + std::to_string(momentumInfNorm_));
    Output::print("Simple", "Maximum cell continuity error : " + std::to_string(maxContinuityError) + "\n");

    momentumBiCGStabIters_= 0;
    pCorrBiCGStabIters_ = 0;
}
