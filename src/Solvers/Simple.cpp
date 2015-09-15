/**
 * @file    Simple.cpp
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
#include "FvScalarScheme.h"
#include "FvVectorScheme.h"
#include "InitialConditions.h"

// ************* Constructors and Destructors *************

Simple::Simple(const Input &input, const HexaFvmMesh &mesh)
    :
      Solver::Solver(input, mesh),
      uField_(mesh, Field<Vector3D>::CONSERVED, "u"),
      pField_(mesh, Field<double>::CONSERVED, "p"),
      rhoField_(mesh, Field<double>::CONSERVED, "rho"),
      muField_(mesh, Field<double>::CONSERVED, "mu"),
      massFlowField_(mesh, Field<double>::CONSERVED, "massFlow"),
      pCorrField_(mesh, Field<double>::CONSERVED, "pCorr"),
      gradUField_(mesh, Field<Tensor3D>::CONSERVED, "gradU"),
      gradPField_(mesh, Field<Vector3D>::CONSERVED, "gradP"),
      gradPCorrField_(mesh, Field<Vector3D>::CONSERVED, "gradPCorr"),
      uField0_(mesh, Field<Vector3D>::AUXILLARY, "uField0"),
      uFieldStar_(mesh, Field<Vector3D>::AUXILLARY, "uFieldStar"),
      dField_(mesh, Field<double>::CONSERVED, "dField"),
      hField_(mesh, Field<Vector3D>::CONSERVED, "hField"),
      gradScalarField_(mesh, Field<Vector3D>::AUXILLARY, "gradScalarField"),
      gradVectorField_(mesh, Field<Tensor3D>::AUXILLARY, "gradVectorField"),
      flowBcs_(input, uField_, pField_, rhoField_, muField_, hField_, dField_, pCorrField_)
{
    InitialConditions initialConditions;

    Solver::createMatrices(2, 3, 7);
    mesh_.addArray3DToTecplotOutput(uField_.name, uField_.cellData());
    mesh_.addArray3DToTecplotOutput(pField_.name, pField_.cellData());
    mesh_.addArray3DToTecplotOutput(massFlowField_.name, massFlowField_.cellData());

    if(Solver::solutionType_ == Solver::STEADY)
        nInnerIters_ = 1;
    else
        nInnerIters_ = input.caseParameters.get<int>("Solver.numberOfInnerIterations");

    omegaMomentum_ = input.caseParameters.get<double>("Solver.relaxationFactorMomentum");
    omegaPCorr_ = input.caseParameters.get<double>("Solver.relaxationFactorPCorr");

    flowBcs_.setParallelBoundaries(mesh.getAdjProcNoPtr());

    setConstantFields(input);
    initialConditions.setInitialConditions(uField_);
    initialConditions.setInitialConditions(pField_);
}

// ************* Public Methods *************

double Simple::solve(double timeStep)
{
    uField0_ = uField_;

    for(int innerIterNo = 0; innerIterNo < nInnerIters_; ++innerIterNo)
    {
        computeMomentum(timeStep);
        computePCorr();
        correct();
    }

    return computeContinuityError();
}

void Simple::displayUpdateMessage()
{
    Output::print("Simple", "completed iteration.");
    Output::print("Simple", "BiCGStab iterations: " + std::to_string(biCGStabIters_));
    Output::print("Simple", "Max continuity error: " + std::to_string(continuityError_));
    Output::printLine();
}

// ************* Protected Methods *************

void Simple::setConstantFields(const Input &input)
{
    rhoField_.setAll(input.caseParameters.get<double>("Solver.rho"));
    muField_.setAll(input.caseParameters.get<double>("Solver.mu"));
}

void Simple::computeMomentum(double timeStep)
{
    using namespace std;

    int cols[7];
    double a0P, a[7];

    uFieldStar_ = uField_;

    //- Assemble the coefficient matrix, to be used for each velocity component
    time_.tic();
    for(int k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!mesh_.iMap.isActive(i, j, k))
                    continue;
                else if(Solver::solutionType_ == Solver::UNSTEADY)
                    a0P = rhoField_(i, j, k)*mesh_.cellVol(i, j, k)/timeStep;
                else
                    a0P = 0.;

                a[1] =  min(massFlowField_.faceE(i, j, k), 0.) - muField_.faceE(i, j, k)*mesh_.dE(i, j, k);
                a[2] = -max(massFlowField_.faceW(i, j, k), 0.) - muField_.faceW(i, j, k)*mesh_.dW(i, j, k);
                a[3] =  min(massFlowField_.faceN(i, j, k), 0.) - muField_.faceN(i, j, k)*mesh_.dN(i, j, k);
                a[4] = -max(massFlowField_.faceS(i, j, k), 0.) - muField_.faceS(i, j, k)*mesh_.dS(i, j, k);
                a[5] =  min(massFlowField_.faceT(i, j, k), 0.) - muField_.faceT(i, j, k)*mesh_.dT(i, j, k);
                a[6] = -max(massFlowField_.faceB(i, j, k), 0.) - muField_.faceB(i, j, k)*mesh_.dB(i, j, k);

                a[0] = (max(massFlowField_.faceE(i, j, k), 0.) + muField_.faceE(i, j, k)*mesh_.dE(i, j, k)
                        - min(massFlowField_.faceW(i, j, k), 0.) + muField_.faceW(i, j, k)*mesh_.dW(i, j, k)
                        + max(massFlowField_.faceN(i, j, k), 0.) + muField_.faceN(i, j, k)*mesh_.dN(i, j, k)
                        - min(massFlowField_.faceS(i, j, k), 0.) + muField_.faceS(i, j, k)*mesh_.dS(i, j, k)
                        + max(massFlowField_.faceT(i, j, k), 0.) + muField_.faceT(i, j, k)*mesh_.dT(i, j, k)
                        - min(massFlowField_.faceB(i, j, k), 0.) + muField_.faceB(i, j, k)*mesh_.dB(i, j, k)
                        + a0P)/omegaMomentum_;

                Vector3D b = a0P*uField0_(i, j, k);
                b += -mesh_.cellVol(i, j, k)*gradPField_(i, j, k);
                b += (1. - omegaMomentum_)*a[0]*uFieldStar_(i, j, k);

                b += muField_.faceE(i, j, k)*dot(gradUField_.faceE(i, j, k), mesh_.cE(i, j, k));
                b += muField_.faceW(i, j, k)*dot(gradUField_.faceW(i, j, k), mesh_.cW(i, j, k));
                b += muField_.faceN(i, j, k)*dot(gradUField_.faceN(i, j, k), mesh_.cN(i, j, k));
                b += muField_.faceS(i, j, k)*dot(gradUField_.faceS(i, j, k), mesh_.cS(i, j, k));
                b += muField_.faceT(i, j, k)*dot(gradUField_.faceT(i, j, k), mesh_.cT(i, j, k));
                b += muField_.faceB(i, j, k)*dot(gradUField_.faceB(i, j, k), mesh_.cB(i, j, k));

                flowBcs_.uFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, b);

                int rowNo = mesh_.iMap(i, j, k, 0);
                cols[0] = rowNo;
                cols[1] = mesh_.iMap(i + 1, j, k, 0);
                cols[2] = mesh_.iMap(i - 1, j, k, 0);
                cols[3] = mesh_.iMap(i, j + 1, k, 0);
                cols[4] = mesh_.iMap(i, j - 1, k, 0);
                cols[5] = mesh_.iMap(i, j, k + 1, 0);
                cols[6] = mesh_.iMap(i, j, k - 1, 0);

                A_[0].setRow(rowNo, 7, cols, a);
                b_[0].setValue(rowNo, b.x);
                b_[1].setValue(rowNo, b.y);
                b_[2].setValue(rowNo, b.z);

                dField_(i, j, k) = mesh_.cellVol(i, j, k)/a[0];
            }
        }
    }

    time_.toc();
    Output::print("Simple", "momentum matrix assembly completed in " + time_.elapsedTime());

    time_.tic();
    for(int componentNo = 0, biCGStabIters_ = 0; componentNo < 3; ++componentNo)
    {
        biCGStabIters_ += A_[0].solve(b_[componentNo], x_[componentNo]);

        for(int k = 0; k < mesh_.nCellsK(); ++k)
        {
            for(int j = 0; j < mesh_.nCellsJ(); ++j)
            {
                for(int i = 0; i < mesh_.nCellsI(); ++i)
                {
                    if(mesh_.iMap.isInactive(i, j, k))
                        continue;

                    uField_(i, j, k)(componentNo) = x_[componentNo](mesh_.iMap(i, j, k, 0));

                    if(mesh_.iMap.isActive(i, j, k))
                        hField_(i, j, k)(componentNo) = uField_(i, j, k)(componentNo) + dField_(i, j, k)*gradPField_(i, j, k)(componentNo);
                    else
                        hField_(i, j, k)(componentNo) = uField_(i, j, k)(componentNo);
                }
            }
        }
    }

    time_.toc();
    Output::print("Simple", "momentum matrix solution completed in " + time_.elapsedTime());
    Output::print("Simple", "biCGStab iterations: " + to_string(biCGStabIters_));

    flowBcs_.uFieldBcs.setBoundaries();
    flowBcs_.dFieldBcs.setBoundaries();
    flowBcs_.hFieldBcs.setBoundaries();

    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, dField_);
    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, hField_);

    rhieChowInterpolateFaces();
}

void Simple::computePCorr()
{
    using namespace std;

    int cols[7];
    double a[7];

    //- Assemble the coefficient matrix
    time_.tic();
    for(int k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!mesh_.iMap.isActive(i, j, k))
                    continue;

                a[1] = rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*mesh_.dE(i, j, k);
                a[2] = rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*mesh_.dW(i, j, k);
                a[3] = rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*mesh_.dN(i, j, k);
                a[4] = rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*mesh_.dS(i, j, k);
                a[5] = rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*mesh_.dT(i, j, k);
                a[6] = rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*mesh_.dB(i, j, k);

                a[0] = -(a[1] + a[2] + a[3] + a[4] + a[5] + a[6]);

                double b = massFlowField_.faceE(i, j, k) - massFlowField_.faceW(i, j, k)
                        + massFlowField_.faceN(i, j, k) - massFlowField_.faceS(i, j, k)
                        + massFlowField_.faceT(i, j, k) - massFlowField_.faceB(i, j, k);

                flowBcs_.pCorrFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, b);

                int rowNo = mesh_.iMap(i, j, k, 0);
                cols[0] = mesh_.iMap(i, j, k, 0);
                cols[1] = mesh_.iMap(i + 1, j, k, 0);
                cols[2] = mesh_.iMap(i - 1, j, k, 0);
                cols[3] = mesh_.iMap(i, j + 1, k, 0);
                cols[4] = mesh_.iMap(i, j - 1, k, 0);
                cols[5] = mesh_.iMap(i, j, k + 1, 0);
                cols[6] = mesh_.iMap(i, j, k - 1, 0);

                A_[1].setRow(rowNo, 7, cols, a);
                b_[0].setValue(rowNo, b);
            }
        }
    }

    time_.toc();
    Output::print("Simple", "pressure correction matrix assembly completed in " + time_.elapsedTime());

    time_.tic();
    biCGStabIters_ = A_[1].solve(b_[0], x_[0]);

    for(int k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(mesh_.iMap.isInactive(i, j,k))
                    continue;

                pCorrField_(i, j, k) = x_[0](mesh_.iMap(i, j, k, 0));
            }
        }
    }

    time_.toc();
    Output::print("Simple", "pressure correction matrix solution completed in " + time_.elapsedTime());
    Output::print("Simple", "biCGStab iterations: " + to_string(biCGStabIters_));

    flowBcs_.pCorrFieldBcs.setBoundaries();

    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, pCorrField_);
    FvScalarScheme::computeCellCenteredGradients(FvScalarScheme::DIVERGENCE_THEOREM, pCorrField_, gradPCorrField_);
}

void Simple::correct()
{
    for(int k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(mesh_.iMap.isInactive(i, j, k))
                    continue;

                if(i == 0 && flowBcs_.massCorrectionRequiredWest())
                {
                    massFlowField_.faceW(i, j, k) += -rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*(pCorrField_(i - 1, j, k) - pCorrField_(i, j, k))*mesh_.dW(i, j, k);
                }
                if(j == 0 && flowBcs_.massCorrectionRequiredSouth())
                {
                    massFlowField_.faceS(i, j, k) += -rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*(pCorrField_(i, j - 1, k) - pCorrField_(i, j, k))*mesh_.dS(i, j, k);
                }
                if(k == 0 && flowBcs_.massCorrectionRequiredBottom())
                {
                    massFlowField_.faceB(i, j, k) += -rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*(pCorrField_(i, j, k - 1) - pCorrField_(i, j, k))*mesh_.dB(i, j, k);
                }

                if(i < mesh_.uCellI() || flowBcs_.massCorrectionRequiredEast())
                {
                    massFlowField_.faceE(i, j, k) += -rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*(pCorrField_(i + 1, j, k) - pCorrField_(i, j, k))*mesh_.dE(i, j, k);
                }
                if(j < mesh_.uCellJ() || flowBcs_.massCorrectionRequiredNorth())
                {
                    massFlowField_.faceN(i, j, k) += -rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*(pCorrField_(i, j + 1, k) - pCorrField_(i, j, k))*mesh_.dN(i, j, k);
                }
                if(k < mesh_.uCellK() || flowBcs_.massCorrectionRequiredTop())
                {
                    massFlowField_.faceT(i, j, k) += -rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*(pCorrField_(i, j, k + 1) - pCorrField_(i, j, k))*mesh_.dT(i, j, k);
                }

                massFlowField_(i, j, k) = massFlowField_.faceE(i, j, k) - massFlowField_.faceW(i, j, k)
                        + massFlowField_.faceN(i, j, k) - massFlowField_.faceS(i, j, k)
                        + massFlowField_.faceT(i, j, k) - massFlowField_.faceB(i, j, k);

                pField_(i, j, k) += omegaPCorr_*pCorrField_(i, j, k);
                uField_(i, j, k) += -dField_(i, j, k)*gradPCorrField_(i, j, k);
            }
        }
    }

    flowBcs_.pFieldBcs.setBoundaries();
    flowBcs_.uFieldBcs.setBoundaries();

    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, pField_);
    FvVectorScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, uField_);

    FvScalarScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, pField_, gradPField_);
    FvVectorScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, uField_, gradUField_);
}

double Simple::computeContinuityError()
{
    continuityError_ = 0.;
    for(int k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!mesh_.iMap.isActive(i, j, k))
                    continue;

                if(fabs(massFlowField_(i, j, k)) > continuityError_)
                {
                    continuityError_ = fabs(massFlowField_(i, j, k))/mesh_.cellVol(i, j, k);
                }
            }
        }
    }
    continuityError_ = Parallel::max(continuityError_);

    return continuityError_;
}

void Simple::rhieChowInterpolateFaces()
{
    for(int k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(mesh_.iMap.isInactive(i, j, k))
                    continue;

                if(i == 0)
                {
                    massFlowField_.faceW(i, j, k) = -dot(hField_.faceW(i, j, k), mesh_.fAreaNormW(i, j, k)) + dField_.faceW(i, j, k)*(pField_(i - 1, j, k) - pField_(i, j, k))*mesh_.dW(i, j, k);
                }
                if(j == 0)
                {
                    massFlowField_.faceS(i, j, k) = -dot(hField_.faceS(i, j, k), mesh_.fAreaNormS(i, j, k)) + dField_.faceS(i, j, k)*(pField_(i, j - 1, k) - pField_(i, j, k))*mesh_.dS(i, j, k);
                }
                if(k == 0)
                {
                    massFlowField_.faceB(i, j, k) = -dot(hField_.faceB(i, j, k), mesh_.fAreaNormB(i, j, k)) + dField_.faceB(i, j, k)*(pField_(i, j, k - 1) - pField_(i, j, k))*mesh_.dB(i, j, k);
                }

                massFlowField_.faceE(i, j, k) = dot(hField_.faceE(i, j, k), mesh_.fAreaNormE(i, j, k)) - dField_.faceE(i, j, k)*(pField_(i + 1, j, k) - pField_(i, j, k))*mesh_.dE(i, j, k);
                massFlowField_.faceN(i, j, k) = dot(hField_.faceN(i, j, k), mesh_.fAreaNormN(i, j, k)) - dField_.faceN(i, j, k)*(pField_(i, j + 1, k) - pField_(i, j, k))*mesh_.dN(i, j, k);
                massFlowField_.faceT(i, j, k) = dot(hField_.faceT(i, j, k), mesh_.fAreaNormT(i, j, k)) - dField_.faceT(i, j, k)*(pField_(i, j, k + 1) - pField_(i, j, k))*mesh_.dT(i, j, k);
            }
        }
    }
}
