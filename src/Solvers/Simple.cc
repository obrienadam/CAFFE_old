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
#include "FvScalarScheme.h"
#include "FvVectorScheme.h"

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
      gradVectorField_(mesh, Field<Tensor3D>::AUXILLARY, "gradVectorField")
{
    Solver::createMatrices(2, 3, 7);
    mesh_.addArray3DToTecplotOutput(uField_.name, uField_.cellData());
    mesh_.addArray3DToTecplotOutput(pField_.name, pField_.cellData());
    mesh_.addArray3DToTecplotOutput(massFlowField_.name, massFlowField_.cellData());

    if(Solver::solutionType_ == Solver::STEADY)
        nInnerIters_ = 1;
    else
        nInnerIters_ = input.getInt("numberOfInnerIterations");

    omegaMomentum_ = input.getDouble("relaxationFactorMomentum");
    omegaPCorr_ = input.getDouble("relaxationFactorPCorr");

    setBoundaries(input);
    setConstantFields(input);
}

// ************* Public Methods *************

double Simple::solve(double timeStep)
{
    uField0_ = uField_;

    for(int i = 0; i < nInnerIters_; ++i)
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

void Simple::setBoundaries(const Input &input)
{
    //- East bcs
    if(input.getString("boundaryTypeEast") == "inlet" || input.getString("boundaryTypeEast") == "wall")
    {
        uField_.setEastBoundary("fixed", std::stov(input.getString("boundaryRefVectorEast")));
        pField_.setEastBoundary("zeroGradient", 0.);
        pCorrField_.setEastBoundary("zeroGradient", 0.);
        hField_.setEastBoundary("fixed", std::stov(input.getString("boundaryRefVectorEast")));
    }
    else if(input.getString("boundaryTypeEast") == "outlet")
    {
        uField_.setEastBoundary("zeroGradient", Vector3D(0., 0., 0.));
        pField_.setEastBoundary("fixed", 0.);
        pCorrField_.setEastBoundary("fixed", 0.);
        hField_.setEastBoundary("zeroGradient", Vector3D(0., 0., 0.));
    }
    else
    {
        uField_.setEastBoundary(input.getString("boundaryTypeEast"), std::stov(input.getString("boundaryRefVectorEast")));
        pField_.setEastBoundary(input.getString("boundaryTypeEast"), input.getDouble("boundaryRefValueEast"));
        pCorrField_.setEastBoundary(input.getString("boundaryTypeEast"), input.getDouble("boundaryRefValueEast"));
        hField_.setEastBoundary(input.getString("boundaryTypeEast"), std::stov(input.getString("boundaryRefVectorEast")));
    }

    //- West bcs
    if(input.getString("boundaryTypeWest") == "inlet" || input.getString("boundaryTypeWest") == "wall")
    {
        uField_.setWestBoundary("fixed", std::stov(input.getString("boundaryRefVectorWest")));
        pField_.setWestBoundary("zeroGradient", 0.);
        pCorrField_.setWestBoundary("zeroGradient", 0.);
        hField_.setWestBoundary("fixed", std::stov(input.getString("boundaryRefVectorWest")));
    }
    else if(input.getString("boundaryTypeWest") == "outlet")
    {
        uField_.setWestBoundary("zeroGradient", Vector3D(0., 0., 0.));
        pField_.setWestBoundary("fixed", 0.);
        pCorrField_.setWestBoundary("fixed", 0.);
        hField_.setWestBoundary("zeroGradient", Vector3D(0., 0., 0.));
    }
    else
    {
        uField_.setWestBoundary(input.getString("boundaryTypeWest"), std::stov(input.getString("boundaryRefVectorWest")));
        pField_.setWestBoundary(input.getString("boundaryTypeWest"), input.getDouble("boundaryRefValueWest"));
        pCorrField_.setWestBoundary(input.getString("boundaryTypeWest"), input.getDouble("boundaryRefValueWest"));
        hField_.setWestBoundary(input.getString("boundaryTypeWest"), std::stov(input.getString("boundaryRefVectorWest")));
    }

    //- North bcs
    if(input.getString("boundaryTypeNorth") == "inlet" || input.getString("boundaryTypeNorth") == "wall")
    {
        uField_.setNorthBoundary("fixed", std::stov(input.getString("boundaryRefVectorNorth")));
        pField_.setNorthBoundary("zeroGradient", 0.);
        pCorrField_.setNorthBoundary("zeroGradient", 0.);
        hField_.setNorthBoundary("fixed", std::stov(input.getString("boundaryRefVectorNorth")));
    }
    else if(input.getString("boundaryTypeNorth") == "outlet")
    {
        uField_.setNorthBoundary("zeroGradient", Vector3D(0., 0., 0.));
        pField_.setNorthBoundary("fixed", 0.);
        pCorrField_.setNorthBoundary("fixed", 0.);
        hField_.setNorthBoundary("zeroGradient", Vector3D(0., 0., 0.));
    }
    else
    {
        uField_.setNorthBoundary(input.getString("boundaryTypeNorth"), std::stov(input.getString("boundaryRefVectorNorth")));
        pField_.setNorthBoundary(input.getString("boundaryTypeNorth"), input.getDouble("boundaryRefValueNorth"));
        pCorrField_.setNorthBoundary(input.getString("boundaryTypeNorth"), input.getDouble("boundaryRefValueNorth"));
        hField_.setNorthBoundary(input.getString("boundaryTypeNorth"), std::stov(input.getString("boundaryRefVectorNorth")));
    }

    //- South bcs
    if(input.getString("boundaryTypeSouth") == "inlet" || input.getString("boundaryTypeSouth") == "wall")
    {
        uField_.setSouthBoundary("fixed", std::stov(input.getString("boundaryRefVectorSouth")));
        pField_.setSouthBoundary("zeroGradient", 0.);
        pCorrField_.setSouthBoundary("zeroGradient", 0.);
        hField_.setSouthBoundary("fixed", std::stov(input.getString("boundaryRefVectorSouth")));
    }
    else if(input.getString("boundaryTypeSouth") == "outlet")
    {
        uField_.setSouthBoundary("zeroGradient", Vector3D(0., 0., 0.));
        pField_.setSouthBoundary("fixed", 0.);
        pCorrField_.setSouthBoundary("fixed", 0.);
        hField_.setSouthBoundary("zeroGradient", Vector3D(0., 0., 0.));
    }
    else
    {
        uField_.setSouthBoundary(input.getString("boundaryTypeSouth"), std::stov(input.getString("boundaryRefVectorSouth")));
        pField_.setSouthBoundary(input.getString("boundaryTypeSouth"), input.getDouble("boundaryRefValueSouth"));
        pCorrField_.setSouthBoundary(input.getString("boundaryTypeSouth"), input.getDouble("boundaryRefValueSouth"));
        hField_.setSouthBoundary(input.getString("boundaryTypeSouth"), std::stov(input.getString("boundaryRefVectorSouth")));
    }

    //- Top bcs
    if(input.getString("boundaryTypeTop") == "inlet" || input.getString("boundaryTypeTop") == "wall")
    {
        uField_.setTopBoundary("fixed", std::stov(input.getString("boundaryRefVectorTop")));
        pField_.setTopBoundary("zeroGradient", 0.);
        pCorrField_.setTopBoundary("zeroGradient", 0.);
        hField_.setTopBoundary("fixed", std::stov(input.getString("boundaryRefVectorTop")));
    }
    else if(input.getString("boundaryTypeTop") == "outlet")
    {
        uField_.setTopBoundary("zeroGradient", Vector3D(0., 0., 0.));
        pField_.setTopBoundary("fixed", 0.);
        pCorrField_.setTopBoundary("fixed", 0.);
        hField_.setTopBoundary("zeroGradient", Vector3D(0., 0., 0.));
    }
    else
    {
        uField_.setTopBoundary(input.getString("boundaryTypeTop"), std::stov(input.getString("boundaryRefVectorTop")));
        pField_.setTopBoundary(input.getString("boundaryTypeTop"), input.getDouble("boundaryRefValueTop"));
        pCorrField_.setTopBoundary(input.getString("boundaryTypeTop"), input.getDouble("boundaryRefValueTop"));
        hField_.setTopBoundary(input.getString("boundaryTypeTop"), std::stov(input.getString("boundaryRefVectorTop")));
    }

    //- Bottom bcs
    if(input.getString("boundaryTypeBottom") == "inlet" || input.getString("boundaryTypeBottom") == "wall")
    {
        uField_.setBottomBoundary("fixed", std::stov(input.getString("boundaryRefVectorBottom")));
        pField_.setBottomBoundary("zeroGradient", 0.);
        pCorrField_.setBottomBoundary("zeroGradient", 0.);
        hField_.setBottomBoundary("fixed", std::stov(input.getString("boundaryRefVectorBottom")));
    }
    else if(input.getString("boundaryTypeBottom") == "outlet")
    {
        uField_.setBottomBoundary("zeroGradient", Vector3D(0., 0., 0.));
        pField_.setBottomBoundary("fixed", 0.);
        pCorrField_.setBottomBoundary("fixed", 0.);
        hField_.setBottomBoundary("zeroGradient", Vector3D(0., 0., 0.));
    }
    else
    {
        uField_.setBottomBoundary(input.getString("boundaryTypeBottom"), std::stov(input.getString("boundaryRefVectorBottom")));
        pField_.setBottomBoundary(input.getString("boundaryTypeBottom"), input.getDouble("boundaryRefValueBottom"));
        pCorrField_.setBottomBoundary(input.getString("boundaryTypeBottom"), input.getDouble("boundaryRefValueBottom"));
        hField_.setBottomBoundary(input.getString("boundaryTypeBottom"), std::stov(input.getString("boundaryRefVectorTop")));
    }
}

void Simple::setConstantFields(const Input &input)
{
    int i, j, k;

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                rhoField_(i, j, k) = input.getDouble("rho");
                muField_(i, j, k) = input.getDouble("mu");
            }
        }
    }

    rhoField_.setEastBoundary("zeroGradient", 0.);
    rhoField_.setWestBoundary("zeroGradient", 0.);
    rhoField_.setNorthBoundary("zeroGradient", 0.);
    rhoField_.setSouthBoundary("zeroGradient", 0.);
    rhoField_.setTopBoundary("zeroGradient", 0.);
    rhoField_.setBottomBoundary("zeroGradient", 0.);
    muField_.setEastBoundary("zeroGradient", 0.);
    muField_.setWestBoundary("zeroGradient", 0.);
    muField_.setNorthBoundary("zeroGradient", 0.);
    muField_.setSouthBoundary("zeroGradient", 0.);
    muField_.setTopBoundary("zeroGradient", 0.);
    muField_.setBottomBoundary("zeroGradient", 0.);
    dField_.setEastBoundary("zeroGradient", 0.);
    dField_.setWestBoundary("zeroGradient", 0.);
    dField_.setNorthBoundary("zeroGradient", 0.);
    dField_.setSouthBoundary("zeroGradient", 0.);
    dField_.setTopBoundary("zeroGradient", 0.);
    dField_.setBottomBoundary("zeroGradient", 0.);
    Output::printLine();

    rhoField_.setBoundaryFields();
    muField_.setBoundaryFields();

    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, rhoField_);
    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, muField_);
}

void Simple::computeMomentum(double timeStep)
{
    using namespace std;

    int i, j, k, componentNo, cols[7], rowNo;
    double a0P, a[7];
    Vector3D b = Vector3D(0., 0., 0.);

    uFieldStar_ = uField_;

    //- Assemble the coefficient matrix, to be used for each velocity component
    time_.tic();
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!indexMap_.isActive(i, j, k))
                    continue;
                else if(Solver::solutionType_ == Solver::UNSTEADY)
                    a0P = rhoField_(i, j, k)*mesh_.cellVol(i, j, k)/timeStep;
                else
                    a0P = 0.;

                a[1] =  min(massFlowField_.faceE(i, j, k), 0.) - muField_.faceE(i, j, k)*dE_(i, j, k);
                a[2] = -max(massFlowField_.faceW(i, j, k), 0.) - muField_.faceW(i, j, k)*dW_(i, j, k);
                a[3] =  min(massFlowField_.faceN(i, j, k), 0.) - muField_.faceN(i, j, k)*dN_(i, j, k);
                a[4] = -max(massFlowField_.faceS(i, j, k), 0.) - muField_.faceS(i, j, k)*dS_(i, j, k);
                a[5] =  min(massFlowField_.faceT(i, j, k), 0.) - muField_.faceT(i, j, k)*dT_(i, j, k);
                a[6] = -max(massFlowField_.faceB(i, j, k), 0.) - muField_.faceB(i, j, k)*dB_(i, j, k);

                a[0] = (max(massFlowField_.faceE(i, j, k), 0.) + muField_.faceE(i, j, k)*dE_(i, j, k)
                        - min(massFlowField_.faceW(i, j, k), 0.) + muField_.faceW(i, j, k)*dW_(i, j, k)
                        + max(massFlowField_.faceN(i, j, k), 0.) + muField_.faceN(i, j, k)*dN_(i, j, k)
                        - min(massFlowField_.faceS(i, j, k), 0.) + muField_.faceS(i, j, k)*dS_(i, j, k)
                        + max(massFlowField_.faceT(i, j, k), 0.) + muField_.faceT(i, j, k)*dT_(i, j, k)
                        - min(massFlowField_.faceB(i, j, k), 0.) + muField_.faceB(i, j, k)*dB_(i, j, k)
                        + a0P)/omegaMomentum_;

                b = a0P*uField0_(i, j, k);
                uField_.setImplicitBoundaryCoeffs(i, j, k, a, b);
                b += -mesh_.cellVol(i, j, k)*gradPField_(i, j, k);
                b += (1. - omegaMomentum_)*a[0]*uFieldStar_(i, j, k);

                b += muField_.faceE(i, j, k)*dot(gradUField_.faceE(i, j, k), cE_(i, j, k));
                b += muField_.faceW(i, j, k)*dot(gradUField_.faceW(i, j, k), cW_(i, j, k));
                b += muField_.faceN(i, j, k)*dot(gradUField_.faceN(i, j, k), cN_(i, j, k));
                b += muField_.faceS(i, j, k)*dot(gradUField_.faceS(i, j, k), cS_(i, j, k));
                b += muField_.faceT(i, j, k)*dot(gradUField_.faceT(i, j, k), cT_(i, j, k));
                b += muField_.faceB(i, j, k)*dot(gradUField_.faceB(i, j, k), cB_(i, j, k));

                rowNo = indexMap_(i, j, k, 0);
                cols[0] = rowNo;
                cols[1] = indexMap_(i + 1, j, k, 0);
                cols[2] = indexMap_(i - 1, j, k, 0);
                cols[3] = indexMap_(i, j + 1, k, 0);
                cols[4] = indexMap_(i, j - 1, k, 0);
                cols[5] = indexMap_(i, j, k + 1, 0);
                cols[6] = indexMap_(i, j, k - 1, 0);

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

    biCGStabIters_ = 0;
    time_.tic();
    for(componentNo = 0; componentNo < 3; ++componentNo)
    {
        biCGStabIters_ += A_[0].solve(b_[componentNo], x_[componentNo]);

        for(k = 0; k < mesh_.nCellsK(); ++k)
        {
            for(j = 0; j < mesh_.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh_.nCellsI(); ++i)
                {
                    if(indexMap_.isInactive(i, j, k))
                        continue;

                    uField_(i, j, k)(componentNo) = x_[componentNo](indexMap_(i, j, k, 0));

                    if(indexMap_.isActive(i, j, k))
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

    dField_.setBoundaryFields();
    hField_.setBoundaryFields();

    FvScalarScheme::extrapolateInteriorFaces(FvScheme::DIVERGENCE_THEOREM, dField_, gradScalarField_);
    FvVectorScheme::extrapolateInteriorFaces(FvScheme::DIVERGENCE_THEOREM, hField_, gradVectorField_);

    uField_.setBoundaryFields();
    rhieChowInterpolateFaces();
}

void Simple::computePCorr()
{
    using namespace std;

    int i, j, k, cols[7], rowNo;
    double a[7];
    double b = 0.;

    //- Assemble the coefficient matrix
    time_.tic();
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!indexMap_.isActive(i, j, k))
                    continue;

                a[1] = rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*dE_(i, j, k);
                a[2] = rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*dW_(i, j, k);
                a[3] = rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*dN_(i, j, k);
                a[4] = rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*dS_(i, j, k);
                a[5] = rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*dT_(i, j, k);
                a[6] = rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*dB_(i, j, k);

                a[0] = -(a[1] + a[2] + a[3] + a[4] + a[5] + a[6]);

                b = massFlowField_.faceE(i, j, k) - massFlowField_.faceW(i, j, k)
                        + massFlowField_.faceN(i, j, k) - massFlowField_.faceS(i, j, k)
                        + massFlowField_.faceT(i, j, k) - massFlowField_.faceB(i, j, k);

                pCorrField_.setImplicitBoundaryCoeffs(i, j, k, a, b);

                rowNo = indexMap_(i, j, k, 0);
                cols[0] = indexMap_(i, j, k, 0);
                cols[1] = indexMap_(i + 1, j, k, 0);
                cols[2] = indexMap_(i - 1, j, k, 0);
                cols[3] = indexMap_(i, j + 1, k, 0);
                cols[4] = indexMap_(i, j - 1, k, 0);
                cols[5] = indexMap_(i, j, k + 1, 0);
                cols[6] = indexMap_(i, j, k - 1, 0);

                A_[1].setRow(rowNo, 7, cols, a);
                b_[0].setValue(rowNo, b);
            }
        }
    }

    time_.toc();
    Output::print("Simple", "pressure correction matrix assembly completed in " + time_.elapsedTime());

    time_.tic();
    biCGStabIters_ = A_[1].solve(b_[0], x_[0]);

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(indexMap_.isInactive(i, j,k))
                    continue;

                pCorrField_(i, j, k) = x_[0](indexMap_(i, j, k, 0));
            }
        }
    }

    time_.toc();
    Output::print("Simple", "pressure correction matrix solution completed in " + time_.elapsedTime());
    Output::print("Simple", "biCGStab iterations: " + to_string(biCGStabIters_));

    pCorrField_.setBoundaryFields();
    FvScalarScheme::extrapolateInteriorFaces(FvScheme::DIVERGENCE_THEOREM, pCorrField_, gradPCorrField_);
    FvScalarScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, pCorrField_, gradPCorrField_);
}

void Simple::correct()
{
    int i, j, k;

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(indexMap_.isInactive(i, j, k))
                    continue;

                if(i == 0 && uField_.getWestBoundaryPatch() == Field<Vector3D>::ZERO_GRADIENT)
                {
                    massFlowField_.faceW(i, j, k) += -rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*(pCorrField_(i - 1, j, k) - pCorrField_(i, j, k))*dW_(i, j, k);
                }
                if(j == 0 && uField_.getSouthBoundaryPatch() == Field<Vector3D>::ZERO_GRADIENT)
                {
                    massFlowField_.faceS(i, j, k) += -rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*(pCorrField_(i, j - 1, k) - pCorrField_(i, j, k))*dS_(i, j, k);
                }
                if(k == 0 && uField_.getBottomBoundaryPatch() == Field<Vector3D>::ZERO_GRADIENT)
                {
                    massFlowField_.faceB(i, j, k) += -rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*(pCorrField_(i, j, k - 1) - pCorrField_(i, j, k))*dB_(i, j, k);
                }

                if(i < mesh_.uCellI() || uField_.getEastBoundaryPatch() == Field<Vector3D>::ZERO_GRADIENT)
                {
                    massFlowField_.faceE(i, j, k) += -rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*(pCorrField_(i + 1, j, k) - pCorrField_(i, j, k))*dE_(i, j, k);
                }
                if(j < mesh_.uCellJ() || uField_.getNorthBoundaryPatch() == Field<Vector3D>::ZERO_GRADIENT)
                {
                    massFlowField_.faceN(i, j, k) += -rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*(pCorrField_(i, j + 1, k) - pCorrField_(i, j, k))*dN_(i, j, k);
                }
                if(k < mesh_.uCellK() || uField_.getTopBoundaryPatch() == Field<Vector3D>::ZERO_GRADIENT)
                {
                    massFlowField_.faceT(i, j, k) += -rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*(pCorrField_(i, j, k + 1) - pCorrField_(i, j, k))*dT_(i, j, k);
                }

                massFlowField_(i, j, k) = massFlowField_.faceE(i, j, k) - massFlowField_.faceW(i, j, k)
                        + massFlowField_.faceN(i, j, k) - massFlowField_.faceS(i, j, k)
                        + massFlowField_.faceT(i, j, k) - massFlowField_.faceB(i, j, k);

                pField_(i, j, k) += omegaPCorr_*pCorrField_(i, j, k);
                uField_(i, j, k) += -dField_(i, j, k)*gradPCorrField_(i, j, k);
            }
        }
    }

    pField_.setBoundaryFields();
    uField_.setBoundaryFields();

    FvScalarScheme::extrapolateInteriorFaces(FvScheme::DIVERGENCE_THEOREM, pField_, gradPField_);
    FvVectorScheme::extrapolateInteriorFaces(FvScheme::DIVERGENCE_THEOREM, uField_, gradUField_);

    FvScalarScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, pField_, gradPField_);
    FvVectorScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, uField_, gradUField_);
}

double Simple::computeContinuityError()
{
    int i, j, k;

    continuityError_ = 0.;
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!indexMap_.isActive(i, j, k))
                    continue;

                if(fabs(massFlowField_(i, j, k)) > continuityError_)
                    continuityError_ = fabs(massFlowField_(i, j, k))/mesh_.cellVol(i, j, k);
            }
        }
    }

    return continuityError_;
}

void Simple::rhieChowInterpolateFaces()
{
    int i, j, k;

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(indexMap_.isInactive(i, j, k))
                    continue;

                if(i == 0)
                {
                    massFlowField_.faceW(i, j, k) = -dot(hField_.faceW(i, j, k), mesh_.fAreaNormW(i, j, k)) - dField_.faceW(i, j, k)*(pField_(i - 1, j, k) - pField_(i, j, k))*dW_(i, j, k);
                }
                if(j == 0)
                {
                    massFlowField_.faceS(i, j, k) = -dot(hField_.faceS(i, j, k), mesh_.fAreaNormS(i, j, k)) - dField_.faceS(i, j, k)*(pField_(i, j - 1, k) - pField_(i, j, k))*dS_(i, j, k);
                }
                if(k == 0)
                {
                    massFlowField_.faceB(i, j, k) = -dot(hField_.faceB(i, j, k), mesh_.fAreaNormB(i, j, k)) - dField_.faceB(i, j, k)*(pField_(i, j, k - 1) - pField_(i, j, k))*dB_(i, j, k);
                }

                massFlowField_.faceE(i, j, k) = dot(hField_.faceE(i, j, k), mesh_.fAreaNormE(i, j, k)) - dField_.faceE(i, j, k)*(pField_(i + 1, j, k) - pField_(i, j, k))*dE_(i, j, k);
                massFlowField_.faceN(i, j, k) = dot(hField_.faceN(i, j, k), mesh_.fAreaNormN(i, j, k)) - dField_.faceN(i, j, k)*(pField_(i, j + 1, k) - pField_(i, j, k))*dN_(i, j, k);
                massFlowField_.faceT(i, j, k) = dot(hField_.faceT(i, j, k), mesh_.fAreaNormT(i, j, k)) - dField_.faceT(i, j, k)*(pField_(i, j, k + 1) - pField_(i, j, k))*dT_(i, j, k);
            }
        }
    }
}
