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

#include "Simple.h"
#include "InputStringProcessing.h"

// ************* Constructors and Destructors *************

Simple::Simple()
    :
      uStar_("uStar", PRIMITIVE),
      a0P_("a0P", AUXILLARY),
      aP_("aP", AUXILLARY),
      aE_("aE", AUXILLARY),
      aW_("aW", AUXILLARY),
      aN_("aN", AUXILLARY),
      aS_("aS", AUXILLARY),
      aT_("aT", AUXILLARY),
      aB_("aB", AUXILLARY),
      bP_("bP", AUXILLARY),
      massFlow_("massFlow", PRIMITIVE),
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
      rho_(998.),
      mu_(0.1),
      momentumSorToler_(0.01),
      pCorrSorToler_(0.01),
      maxPCorrSorIters_(200),
      maxMomentumSorIters_(50),
      sorOmega_(1.91)
{
    nu_ = mu_/rho_;
}

// ************* Private Methods *************

void Simple::storeUField(Field<Vector3D> &uField)
{
    int i, j, k;

    for(k = 0; k < nFacesK_; ++k)
    {
        for(j = 0; j < nFacesJ_; ++j)
        {
            for(i = 0; i < nFacesI_; ++i)
            {
                if(j < uFaceJ_ && k < uFaceK_)
                    uStar_.faceI(i, j, k) = uField.faceI(i, j, k);

                if(i < uFaceI_ && k < uFaceK_)
                    uStar_.faceJ(i, j, k) = uField.faceJ(i, j, k);

                if(i < uFaceI_ && j < uFaceJ_)
                    uStar_.faceK(i, j, k) = uField.faceK(i, j, k);
            }
        }
    }

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                uStar_(i, j, k) = uField(i, j, k);
            }
        }
    }
}

void Simple::computeMomentum(double timeStep, Field<Vector3D>& uField, Field<double>& pField)
{
    using namespace std;

    int i, j, k, l, itrNo;
    HexaFvmMesh& mesh = *meshPtr_;
    double dE, dW, dN, dS, dT, dB;
    Vector3D sf, ds;
    Vector3D old;

    storeUField(uField);

    interpolateInteriorFaces(uField, VOLUME_WEIGHTED);
    interpolateInteriorFaces(pField, VOLUME_WEIGHTED);

    computeCellCenteredGradients(uField, gradUField_, DIVERGENCE_THEOREM);
    computeCellCenteredGradients(pField, gradPField_, DIVERGENCE_THEOREM);

    interpolateInteriorFaces(gradUField_, VOLUME_WEIGHTED);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                // Time term

                if(true)
                    a0P_(i, j, k) = rho_*mesh.cellVol(i, j, k)/timeStep;

                // Diffusion terms (these are technically constants and could be stored)
                dE = mu_*dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                dW = mu_*dot(mesh.fAreaNormW(i, j, k), mesh.fAreaNormW(i, j, k))/dot(mesh.fAreaNormW(i, j, k), mesh.rCellW(i, j, k));
                dN = mu_*dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                dS = mu_*dot(mesh.fAreaNormS(i, j, k), mesh.fAreaNormS(i, j, k))/dot(mesh.fAreaNormS(i, j, k), mesh.rCellS(i, j, k));
                dT = mu_*dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                dB = mu_*dot(mesh.fAreaNormB(i, j, k), mesh.fAreaNormB(i, j, k))/dot(mesh.fAreaNormB(i, j, k), mesh.rCellB(i, j, k));

                // Face coefficients
                aE_(i, j, k) =  min(massFlow_.faceE(i, j, k), 0.) - dE;
                aW_(i, j, k) =  -max(massFlow_.faceW(i, j, k), 0.) - dW;
                aN_(i, j, k) =  min(massFlow_.faceN(i, j, k), 0.) - dN;
                aS_(i, j, k) =  -max(massFlow_.faceS(i, j, k), 0.) - dS;
                aT_(i, j, k) =  min(massFlow_.faceT(i, j, k), 0.) - dT;
                aB_(i, j, k) =  -max(massFlow_.faceB(i, j, k), 0.) - dB;

                // Central coefficient
                aP_(i, j, k) = (max(massFlow_.faceE(i, j, k), 0.) + dE
                                - min(massFlow_.faceW(i, j, k), 0.) + dW
                                + max(massFlow_.faceN(i, j, k), 0.) + dN
                                - min(massFlow_.faceS(i, j, k), 0.) + dS
                                + max(massFlow_.faceT(i, j, k), 0.) + dT
                                - min(massFlow_.faceB(i, j, k), 0.) + dB
                                + a0P_(i, j, k))/relaxationFactorMomentum_;

                bP_(i, j, k) = a0P_(i, j, k)*uStar_(i, j, k);

                // Compute the cross-diffusion terms (due to mesh non-orthogonality)
/*
                getFaceStencil(i, j, k, EAST, sf, ds);
                bP_(i, j, k) += mu_*dot(gradUField_.faceE(i, j, k), sf - ds*dot(sf, sf)/dot(sf, ds));

                getFaceStencil(i, j, k, WEST, sf, ds);
                bP_(i, j, k) += mu_*dot(gradUField_.faceW(i, j, k), sf - ds*dot(sf, sf)/dot(sf, ds));

                getFaceStencil(i, j, k, NORTH, sf, ds);
                bP_(i, j, k) += mu_*dot(gradUField_.faceN(i, j, k), sf - ds*dot(sf, sf)/dot(sf, ds));

                getFaceStencil(i, j, k, SOUTH, sf, ds);
                bP_(i, j, k) += mu_*dot(gradUField_.faceS(i, j, k), sf - ds*dot(sf, sf)/dot(sf, ds));

                getFaceStencil(i, j, k, TOP, sf, ds);
                bP_(i, j, k) += mu_*dot(gradUField_.faceT(i, j, k), sf - ds*dot(sf, sf)/dot(sf, ds));

                getFaceStencil(i, j, k, BOTTOM, sf, ds);
                bP_(i, j, k) += mu_*dot(gradUField_.faceB(i, j, k), sf - ds*dot(sf, sf)/dot(sf, ds));
*/
                // Higher order terms go here

                bP_(i, j, k) += -mesh.cellVol(i, j, k)*gradPField_(i, j, k);

                // Relaxation source term

                bP_(i, j, k) += (1. - relaxationFactorMomentum_)*aP_(i, j, k)*uStar_(i, j, k);

                // Update D-field

                dField_(i, j, k) = mesh.cellVol(i, j, k)/aP_(i, j, k);
            }
        }
    }

    dField_.setBoundaryFields();

    momentumSorItrs_ = 0;

    for(itrNo = 1; itrNo <= maxMomentumSorIters_; ++itrNo)
    {
        momentumSorConvergence_ = 0.;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    old = uField(i, j, k);

                    for(l = 0; l < 3; ++l)
                    {
                        uField(i, j, k)(l) = (1. - sorOmega_)*uField(i, j, k)(l) + sorOmega_*(- aE_(i, j, k)*uField(i + 1, j, k)(l)
                                                                                              - aW_(i, j, k)*uField(i - 1, j, k)(l)
                                                                                              - aN_(i, j, k)*uField(i, j + 1, k)(l)
                                                                                              - aS_(i, j, k)*uField(i, j - 1, k)(l)
                                                                                              - aT_(i, j, k)*uField(i, j, k + 1)(l)
                                                                                              - aB_(i, j, k)*uField(i, j, k - 1)(l)
                                                                                              + bP_(i, j, k)(l))/aP_(i, j, k);
                    }

                    momentumSorConvergence_ = max((uField(i, j, k) - old).mag(), momentumSorConvergence_);
                }
            }
        }

        uField.setBoundaryFields();

        ++momentumSorItrs_;

        if(momentumSorConvergence_ < momentumSorToler_)
            break;
    }

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
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

    momentumResidual_ = computeResidual(uStar_);
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
                if(i < uCellI_)
                {
                    getFaceStencil(i, j, k, EAST, sf, ds);
                    uField.faceE(i, j, k) = hField_.faceE(i, j, k) - dField_.faceE(i, j, k)*(pField(i + 1, j, k) - pField(i, j, k))*sf/dot(sf, ds);
                }

                if(j < uCellJ_)
                {
                    getFaceStencil(i, j, k, NORTH, sf, ds);
                    uField.faceN(i, j, k) = hField_.faceN(i, j, k) - dField_.faceN(i, j, k)*(pField(i, j + 1, k) - pField(i, j, k))*sf/dot(sf, ds);
                }

                if(k < uCellK_)
                {
                    getFaceStencil(i, j, k, TOP, sf, ds);
                    uField.faceT(i, j, k) = hField_.faceT(i, j, k) - dField_.faceT(i, j, k)*(pField(i, j, k + 1) - pField(i, j, k))*sf/dot(sf, ds);
                }
            }
        }
    }
}

void Simple::computeMassFlowFaces(Field<Vector3D> &uField)
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
                    massFlow_.faceI(i, j, k) = rho_*dot(uField.faceI(i, j, k), mesh.fAreaNormI(i, j, k));

                if(i < uFaceI_ && k < uFaceK_)
                    massFlow_.faceJ(i, j, k) = rho_*dot(uField.faceJ(i, j, k), mesh.fAreaNormJ(i, j, k));

                if(i < uFaceI_ && j < uFaceJ_)
                    massFlow_.faceK(i, j, k) = rho_*dot(uField.faceK(i, j, k), mesh.fAreaNormK(i, j, k));
            }
        }
    }
}

void Simple::computePCorr(Field<Vector3D>& uField, Field<double>& pField)
{
    int i, j, k, itrNo;
    HexaFvmMesh& mesh = *meshPtr_;
    Vector3D sf, ds, gradPCorrBar;
    double old = 0.;

    rhieChowInterpolateInteriorFaces(uField, pField);
    computeMassFlowFaces(uField);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                aE_(i, j, k) = rho_*dField_.faceE(i, j, k)*dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                aW_(i, j, k) = rho_*dField_.faceW(i, j, k)*dot(mesh.fAreaNormW(i, j, k), mesh.fAreaNormW(i, j, k))/dot(mesh.fAreaNormW(i, j, k), mesh.rCellW(i, j, k));
                aN_(i, j, k) = rho_*dField_.faceN(i, j, k)*dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                aS_(i, j, k) = rho_*dField_.faceS(i, j, k)*dot(mesh.fAreaNormS(i, j, k), mesh.fAreaNormS(i, j, k))/dot(mesh.fAreaNormS(i, j, k), mesh.rCellS(i, j, k));
                aT_(i, j, k) = rho_*dField_.faceT(i, j, k)*dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                aB_(i, j, k) = rho_*dField_.faceB(i, j, k)*dot(mesh.fAreaNormB(i, j, k), mesh.fAreaNormB(i, j, k))/dot(mesh.fAreaNormB(i, j, k), mesh.rCellB(i, j, k));

                aP_(i, j, k) = -(aE_(i, j, k) + aW_(i, j, k) + aN_(i, j, k) + aS_(i, j, k) + aT_(i, j, k) + aB_(i, j, k));

                massFlow_(i, j, k) = massFlow_.faceW(i, j, k) - massFlow_.faceE(i, j, k)
                        + massFlow_.faceS(i, j, k) - massFlow_.faceN(i, j, k)
                        + massFlow_.faceB(i, j, k) - massFlow_.faceT(i, j, k);
            }
        }
    }

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                //- Mass flow source

                bP_(i, j, k).x = -massFlow_(i, j, k);
            }
        }
    }

    //- iterate using successive over-relaxation
    pCorrSorItrs_ = 0;

    for(itrNo = 1; itrNo <= maxPCorrSorIters_; ++itrNo)
    {
        pCorrSorConvergence_ = 0.;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    old = pCorr_(i, j, k);

                    pCorr_(i, j, k) = (1. - sorOmega_)*pCorr_(i, j, k) + sorOmega_/aP_(i, j, k)*(bP_(i, j, k).x
                                                                                                 - aE_(i, j, k)*pCorr_(i + 1, j, k)
                                                                                                 - aW_(i, j, k)*pCorr_(i - 1, j, k)
                                                                                                 - aN_(i, j, k)*pCorr_(i, j + 1, k)
                                                                                                 - aS_(i, j, k)*pCorr_(i, j - 1, k)
                                                                                                 - aT_(i, j, k)*pCorr_(i, j, k + 1)
                                                                                                 - aB_(i, j, k)*pCorr_(i, j, k - 1));

                    pCorrSorConvergence_ = std::max(fabs(pCorr_(i, j, k) - old), pCorrSorConvergence_);
                }
            }
        }

        pCorr_.setBoundaryFields();
        ++pCorrSorItrs_;

        if(pCorrSorConvergence_ < pCorrSorToler_)
            break;
    }

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                mesh.findScalarField("pCorr")(i, j, k) = pCorr_(i, j, k);
            }
        }
    }
}

void Simple::correctContinuity(Field<Vector3D> &uField, Field<double> &pField)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& massFlow = mesh.findScalarField("massFlow");
    double a;

    // Correct the mass flow field

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(i < uCellI_)
                {
                    a = -rho_*dField_.faceE(i, j, k)*dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                    massFlow_.faceE(i, j, k) += a*(pCorr_(i + 1, j, k) - pCorr_(i, j, k));
                }

                if(j < uCellJ_)
                {
                    a = -rho_*dField_.faceN(i, j, k)*dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                    massFlow_.faceN(i, j, k) += a*(pCorr_(i, j + 1, k) - pCorr_(i, j, k));
                }

                if(k < uCellK_)
                {
                    a = -rho_*dField_.faceT(i, j, k)*dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                    massFlow_.faceT(i, j, k) += a*(pCorr_(i, j, k + 1) - pCorr_(i, j, k));
                }

                massFlow_(i, j, k) = massFlow_.faceE(i, j, k) - massFlow_.faceW(i, j, k)
                        + massFlow_.faceN(i, j, k) - massFlow_.faceS(i, j, k)
                        + massFlow_.faceT(i, j, k) - massFlow_.faceB(i, j, k);

                massFlow(i, j, k) = fabs(massFlow_(i, j, k));
            }
        }
    }

    // Correct the pressure field

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                pField(i, j, k) += relaxationFactorPCorr_*pCorr_(i, j, k);
            }
        }
    }

    pField.setBoundaryFields();

    // Correct the velocity field

    interpolateInteriorFaces(pCorr_, VOLUME_WEIGHTED);
    computeCellCenteredGradients(pCorr_, gradPCorr_, DIVERGENCE_THEOREM);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                uField(i, j, k) += -dField_(i, j, k)*gradPCorr_(i, j, k);
            }
        }
    }

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

    return sqrt(sumMomentumResidualSqr/sumVol);
}

// ************* Public Methods *************

void Simple::initialize(Input &input, HexaFvmMesh &mesh)
{
    FvScheme::initialize(input, mesh, "NA");
    uFieldPtr_ = &mesh.findVectorField("u");
    pFieldPtr_ = &mesh.findScalarField("p");

    if(input.inputStrings["timeAccurate"] == "ON")
        timeAccurate_ = true;
    else if (input.inputStrings["timeAccurate"] == "OFF")
        timeAccurate_ = false;

    relaxationFactorMomentum_ = input.inputDoubles["relaxationFactorMomentum"];
    relaxationFactorPCorr_ = input.inputDoubles["relaxationFactorPCorr"];
    rho_ = input.inputDoubles["rho"];
    mu_ = input.inputDoubles["mu"];
    momentumSorToler_ = input.inputDoubles["momentumSorToler"];
    pCorrSorToler_ = input.inputDoubles["pCorrSorToler"];
    maxMomentumSorIters_ = input.inputInts["maxMomentumSorIters"];
    maxPCorrSorIters_ = input.inputInts["maxPCorrSorIters"];
    sorOmega_ = input.inputDoubles["sorOmega"];
    nu_ = mu_/rho_;

    uStar_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    a0P_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aP_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aE_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aW_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aN_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aS_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aT_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aB_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    bP_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    hField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dField_.allocate(nCellsI_, nCellsJ_, nCellsK_);


    massFlow_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    pCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    gradPCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradUField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradPField_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    gradVecField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradScalarField_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    setBoundaryConditions(input);
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

        uBoundaryRefValues[i] = Vector3D(buffer);

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

void Simple::discretize(std::vector<double> &timeDerivatives_)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;

    computeMomentum(2e-5, uField, pField);
    computePCorr(uField, pField);
    correctContinuity(uField, pField);
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
                maxContinuityError = std::max(fabs(massFlow_(i, j, k)), maxContinuityError);
            }
        }
    }

    Output::print("Simple", "Momentum prediction SOR iterations  : " + std::to_string(momentumSorItrs_));
    Output::print("Simple", "Momentum prediction SOR convergence : " + std::to_string(momentumSorConvergence_));
    Output::print("Simple", "Pressure correction SOR iterations  : " + std::to_string(pCorrSorItrs_));
    Output::print("Simple", "Pressure correction SOR convergence : " + std::to_string(pCorrSorConvergence_) + "\n");
    Output::print("Simple", "U-Momentum residual                 : " + std::to_string(momentumResidual_.x));
    Output::print("Simple", "V-Momentum residual                 : " + std::to_string(momentumResidual_.y));
    Output::print("Simple", "W-Momentum residual                 : " + std::to_string(momentumResidual_.z));
    Output::print("Simple", "Maximum cell continuity error       : " + std::to_string(maxContinuityError) + "\n");
}
