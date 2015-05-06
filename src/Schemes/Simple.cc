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

// ************* Constructors and Destructors *************

Simple::Simple()
    :
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
      dField_("dField", PRIMITIVE),
      gradUField_("gradUField", PRIMITIVE),
      gradPField_("gradPField", PRIMITIVE),
      relaxationFactorMomentum_(0.3),
      relaxationFactorPCorr_(0.1),
      gradReconstructionMethod_(LEAST_SQUARES),
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

void Simple::rhieChowInterpolateFaces(Field<Vector3D> &uField, Field<double>& pField, Field<double>& dField)
{
    int i, j, k, uI, uJ, uK;
    double alpha;
    Vector3D sf, ds;

    uI = nCellsI_ - 1;
    uJ = nCellsJ_ - 1;
    uK = nCellsK_ - 1;

    computeCellCenteredGradients(pField, gradPField_, gradReconstructionMethod_);
    computeFaceCenteredGradients(pField, gradPField_);

    interpolateInteriorFaces(uField, DISTANCE_WEIGHTED);
    interpolateInteriorFaces(gradPField_, VOLUME_WEIGHTED);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(i < uI)
                {
                    getMeshStencil(i, j, k, EAST, sf, ds, alpha);
                    uField.faceE(i, j, k) += -(dField(i, j, k) + dField(i + 1, j, k))*(pField(i + 1, j, k) - pField(i, j, k) - dot(gradPField_.faceE(i, j, k), ds))*sf/dot(sf, ds);
                }

                if(j < uJ)
                {
                    getMeshStencil(i, j, k, NORTH, sf, ds, alpha);
                    uField.faceN(i, j, k) += -(dField(i, j, k) + dField(i, j + 1, k))*(pField(i, j + 1, k) - pField(i, j, k) - dot(gradPField_.faceN(i, j, k), ds))*sf/dot(sf, ds);
                }

                if(k < uK)
                {
                    getMeshStencil(i, j, k, TOP, sf, ds, alpha);
                    uField.faceT(i, j, k) += -(dField(i, j, k) + dField(i, j, k + 1))*(pField(i, j, k + 1) - pField(i, j, k) - dot(gradPField_.faceT(i, j, k), ds))*sf/dot(sf, ds);
                }
            }
        }
    }
}

void Simple::computeMassFlow(Field<Vector3D> &uField)
{
    int i, j, k, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;

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
                    massFlow_.faceI(i, j, k) = rho_*dot(uField.faceI(i, j, k), mesh.fAreaNormI(i, j, k));

                if(i < uI && k < uK)
                    massFlow_.faceJ(i, j, k) = rho_*dot(uField.faceJ(i, j, k), mesh.fAreaNormJ(i, j, k));

                if(i < uI && j < uJ)
                    massFlow_.faceK(i, j, k) = rho_*dot(uField.faceK(i, j, k), mesh.fAreaNormK(i, j, k));
            }
        }
    }
}

void Simple::computeMomentum(Field<Vector3D>& uField, Field<double>& pField)
{
    using namespace std;

    int i, j, k, l, itrNo;
    HexaFvmMesh& mesh = *meshPtr_;
    double alpha, dE, dW, dN, dS, dT, dB;
    Vector3D sf, ds;
    Tensor3D gradUBar;
    Vector3D old;

    uField.setBoundaryFields();
    pField.setBoundaryFields();

    computeCellCenteredGradients(pField, gradPField_, gradReconstructionMethod_);
    computeCellCenteredJacobians(uField, gradUField_, gradReconstructionMethod_);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
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
                aP_(i, j, k) = max(massFlow_.faceE(i, j, k), 0.) + dE
                        - min(massFlow_.faceW(i, j, k), 0.) + dW
                        + max(massFlow_.faceN(i, j, k), 0.) + dN
                        - min(massFlow_.faceS(i, j, k), 0.) + dS
                        + max(massFlow_.faceT(i, j, k), 0.) + dT
                        - min(massFlow_.faceB(i, j, k), 0.) + dB;

                // Compute the cross-diffusion terms (due to mesh non-orthogonality). These still need to be checked for correct signs

                FvScheme::getMeshStencil(i, j, k, EAST, sf, ds, alpha);
                gradUBar = alpha*gradUField_(i, j, k) + (1. - alpha)*gradUField_(i + 1, j, k);

                bP_(i, j, k) = dot(mu_*gradUBar, sf - ds*dot(sf, sf)/dot(sf, ds));

                FvScheme::getMeshStencil(i, j, k, WEST, sf, ds, alpha);
                gradUBar = alpha*gradUField_(i, j, k) + (1. - alpha)*gradUField_(i - 1, j, k);

                bP_(i, j, k) += dot(mu_*gradUBar, sf - ds*dot(sf, sf)/dot(sf, ds));

                FvScheme::getMeshStencil(i, j, k, NORTH, sf, ds, alpha);
                gradUBar = alpha*gradUField_(i, j, k) + (1. - alpha)*gradUField_(i, j + 1, k);

                bP_(i, j, k) += dot(mu_*gradUBar, sf - ds*dot(sf, sf)/dot(sf, ds));

                FvScheme::getMeshStencil(i, j, k, SOUTH, sf, ds, alpha);
                gradUBar = alpha*gradUField_(i, j, k) + (1. - alpha)*gradUField_(i, j - 1, k);

                bP_(i, j, k) += dot(mu_*gradUBar, sf - ds*dot(sf, sf)/dot(sf, ds));

                FvScheme::getMeshStencil(i, j, k, TOP, sf, ds, alpha);
                gradUBar = alpha*gradUField_(i, j, k) + (1. - alpha)*gradUField_(i, j, k + 1);

                bP_(i, j, k) += dot(mu_*gradUBar, sf - ds*dot(sf, sf)/dot(sf, ds));

                FvScheme::getMeshStencil(i, j, k, BOTTOM, sf, ds, alpha);
                gradUBar = alpha*gradUField_(i, j, k) + (1. - alpha)*gradUField_(i, j, k - 1);

                bP_(i, j, k) += dot(mu_*gradUBar, sf - ds*dot(sf, sf)/dot(sf, ds));

                // Higher order terms go here

                bP_(i, j, k) += -gradPField_(i, j, k)*mesh.cellVol(i, j, k);

                // Relaxation source term

                bP_(i, j, k) += (1. - relaxationFactorMomentum_)/relaxationFactorMomentum_*aP_(i, j, k)*uField(i, j, k);

                // Update D-field

                dField_(i, j, k) = mesh.cellVol(i, j, k)/aP_(i, j, k);
            }
        }
    }

    for(itrNo = 1; itrNo <= maxMomentumSorIters_; ++itrNo)
    {
        uField.setBoundaryFields();
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
                        uField(i, j, k)(l) = (1. - sorOmega_)*uField(i, j, k)(l) + relaxationFactorMomentum_*sorOmega_/aP_(i, j, k)*(bP_(i, j, k)(l)
                                                                                                                                     - aE_(i, j, k)*uField(i + 1, j, k)(l)
                                                                                                                                     - aW_(i, j, k)*uField(i - 1, j, k)(l)
                                                                                                                                     - aN_(i, j, k)*uField(i, j + 1, k)(l)
                                                                                                                                     - aS_(i, j, k)*uField(i, j - 1, k)(l)
                                                                                                                                     - aT_(i, j, k)*uField(i, j, k + 1)(l)
                                                                                                                                     - aB_(i, j, k)*uField(i, j, k - 1)(l));
                    }

                    momentumSorConvergence_ = max((uField(i, j, k) - old).mag(), momentumSorConvergence_);
                }
            }
        }

        if(momentumSorConvergence_ < momentumSorToler_)
            break;
    }
}

void Simple::computePCorr(Field<Vector3D>& uField, Field<double>& pField)
{
    int i, j, k, itrNo;
    HexaFvmMesh& mesh = *meshPtr_;
    Vector3D sf, ds, gradPCorrBar;
    double old = 0.;

    rhieChowInterpolateFaces(uField, pField, dField_);
    computeMassFlow(uField);

    computeCellCenteredGradients(pCorr_, gradPCorr_, gradReconstructionMethod_);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                aE_(i, j, k) = rho_*(dField_(i, j, k) + dField_(i + 1, j, k))*dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                aW_(i, j, k) = rho_*(dField_(i, j, k) + dField_(i - 1, j, k))*dot(mesh.fAreaNormW(i, j, k), mesh.fAreaNormW(i, j, k))/dot(mesh.fAreaNormW(i, j, k), mesh.rCellW(i, j, k));
                aN_(i, j, k) = rho_*(dField_(i, j, k) + dField_(i, j + 1, k))*dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                aS_(i, j, k) = rho_*(dField_(i, j, k) + dField_(i, j - 1, k))*dot(mesh.fAreaNormS(i, j, k), mesh.fAreaNormS(i, j, k))/dot(mesh.fAreaNormS(i, j, k), mesh.rCellS(i, j, k));
                aT_(i, j, k) = rho_*(dField_(i, j, k) + dField_(i, j, k + 1))*dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                aB_(i, j, k) = rho_*(dField_(i, j, k) + dField_(i, j, k - 1))*dot(mesh.fAreaNormB(i, j, k), mesh.fAreaNormB(i, j, k))/dot(mesh.fAreaNormB(i, j, k), mesh.rCellB(i, j, k));

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

    for(itrNo = 1; itrNo <= maxPCorrSorIters_; ++itrNo)
    {
        pCorr_.setBoundaryFields();
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

void Simple::correctMassFlow()
{
    int i, j, k, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;
    Field<double>& massFlow = mesh.findScalarField("massFlow");
    double a;

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
                    a = -rho_*(dField_(i, j, k) + dField_(i + 1, j, k))*dot(mesh.fAreaNormE(i, j, k), mesh.fAreaNormE(i, j, k))/dot(mesh.fAreaNormE(i, j, k), mesh.rCellE(i, j, k));
                    massFlow_.faceE(i, j, k) += a*(pCorr_(i + 1, j, k) - pCorr_(i, j, k));
                }

                if(j < uJ)
                {
                    a = -rho_*(dField_(i, j, k) + dField_(i, j + 1, k))*dot(mesh.fAreaNormN(i, j, k), mesh.fAreaNormN(i, j, k))/dot(mesh.fAreaNormN(i, j, k), mesh.rCellN(i, j, k));
                    massFlow_.faceN(i, j, k) += a*(pCorr_(i, j + 1, k) - pCorr_(i, j, k));
                }

                if(k < uK)
                {
                    a = -rho_*(dField_(i, j, k) + dField_(i, j, k + 1))*dot(mesh.fAreaNormT(i, j, k), mesh.fAreaNormT(i, j, k))/dot(mesh.fAreaNormT(i, j, k), mesh.rCellT(i, j, k));
                    massFlow_.faceT(i, j, k) += a*(pCorr_(i, j, k + 1) - pCorr_(i, j, k));
                }

                massFlow_(i, j, k) = massFlow_.faceE(i, j, k) - massFlow_.faceW(i, j, k)
                        + massFlow_.faceN(i, j, k) - massFlow_.faceS(i, j, k)
                        + massFlow_.faceT(i, j, k) - massFlow_.faceB(i, j, k);

                massFlow(i, j, k) = fabs(massFlow_(i, j, k));
            }
        }
    }
}

void Simple::correctPressure(Field<double> &pField)
{
    int i, j, k;

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
}

void Simple::correctVelocity(Field<Vector3D> &uField)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;

    pCorr_.setBoundaryFields();
    interpolateInteriorFaces(pCorr_, DISTANCE_WEIGHTED);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                uField(i, j, k) += -dField_(i, j, k)/mesh.cellVol(i, j, k)*(pCorr_.faceE(i, j, k)*mesh.fAreaNormE(i, j, k)
                                                                            + pCorr_.faceW(i, j, k)*mesh.fAreaNormW(i, j, k)
                                                                            + pCorr_.faceN(i, j, k)*mesh.fAreaNormN(i, j, k)
                                                                            + pCorr_.faceS(i, j, k)*mesh.fAreaNormS(i, j, k)
                                                                            + pCorr_.faceT(i, j, k)*mesh.fAreaNormT(i, j, k)
                                                                            + pCorr_.faceB(i, j, k)*mesh.fAreaNormB(i, j, k));
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

    aP_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aE_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aW_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aN_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aS_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aT_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    aB_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    bP_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    massFlow_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    pCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradPCorr_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    dField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradUField_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    gradPField_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    pCorr_.setAllBoundaries(pFieldPtr_->getEastBoundaryPatch(), 0.,
                            pFieldPtr_->getWestBoundaryPatch(), 0.,
                            pFieldPtr_->getNorthBoundaryPatch(), 0.,
                            pFieldPtr_->getSouthBoundaryPatch(), 0.,
                            pFieldPtr_->getTopBoundaryPatch(), 0.,
                            pFieldPtr_->getBottomBoundaryPatch(), 0.);

    dField_.setAllBoundaries(ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.,
                             ZERO_GRADIENT, 0.);
}

int Simple::nConservedVariables()
{
    return 4*meshPtr_->size();
}

void Simple::discretize(std::vector<double> &timeDerivatives_)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;

    //- Compute the predicted momentum

    computeMomentum(uField, pField);

    //- Compute the pressure corrections

    computePCorr(uField, pField);

    //- Apply the mass flow, pressure and velocity corrections

    correctMassFlow();
    correctPressure(pField);
    correctVelocity(uField);
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

    Output::print("Simple", "Momentum Prediction SOR Convergence: " + std::to_string(momentumSorConvergence_));
    Output::print("Simple", "Pressure Correction SOR Convergence: " + std::to_string(pCorrSorConvergence_));
    Output::print("Simple", "Maximum continuity error: " + std::to_string(maxContinuityError) + "\n");
}
