#include "Coupled.h"
#include "FvScalarScheme.h"
#include "FvVectorScheme.h"

Coupled::Coupled(const Input &input, const HexaFvmMesh &mesh)
    :
      Solver(input, mesh),
      uField_(mesh, Field<Vector3D>::CONSERVED, "u"),
      pField_(mesh, Field<double>::CONSERVED, "p"),
      rhoField_(mesh, Field<double>::CONSERVED, "rho"),
      muField_(mesh, Field<double>::CONSERVED, "mu"),
      gradPField_(mesh, Field<Vector3D>::CONSERVED, "gradPField"),
      dField_(mesh, Field<double>::CONSERVED, "d"),
      hField_(mesh, Field<Vector3D>::CONSERVED, "h"),
      massFlowField_(mesh, Field<double>::CONSERVED, "massFlow"),
      uField0_(mesh, Field<Vector3D>::AUXILLARY, "uField0"),
      uFieldStar_(mesh, Field<Vector3D>::AUXILLARY, "uFieldStar"),
      flowBcs_(input, uField_, pField_, rhoField_, muField_, hField_, dField_)
{
    Solver::createMatrices(4, 1, 1, 28);
    mesh_.addArray3DToTecplotOutput(uField_.name, uField_.cellData());
    mesh_.addArray3DToTecplotOutput(pField_.name, pField_.cellData());
    mesh_.addArray3DToTecplotOutput(massFlowField_.name, massFlowField_.cellData());

    if(Solver::solutionType_ == Solver::STEADY)
        nInnerIters_ = 1;
    else
        nInnerIters_ = input.caseParameters.get<int>("Solver.numberOfInnerIterations");

    omegaMomentum_ = input.caseParameters.get<double>("Solver.relaxationFactorMomentum");

    rhoField_.setAll(input.caseParameters.get<double>("Solver.rho"));
    muField_.setAll(input.caseParameters.get<double>("Solver.mu"));
}

//*************************** Public methods ********************************

double Coupled::solve(double timeStep)
{
    int innerIterNo;

    uField0_ = uField_;

    for(innerIterNo = 0; innerIterNo < nInnerIters_; ++innerIterNo)
    {
        computeMomentumAndPressure(timeStep);
    }
}

//*************************** Protected methods *****************************

void Coupled::computeMomentumAndPressure(double timeStep)
{
    using namespace std;

    int i, j, k, componentNo, cols[14], rowNo;
    double a0P, a[14], bp = 0.;
    Vector3D bu = Vector3D(0., 0., 0.);

    uFieldStar_ = uField_;

    time_.tic();
    /** Momentum equation **/
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                if(!indexMap_.isActive(i, j, k))
                    continue;
                else if(solutionType_ == UNSTEADY)
                    a0P = rhoField_(i, j, k)*mesh_.cellVol(i, j, k)/timeStep;
                else
                    a0P = 0.;

                //- Compute the advection/diffusion coefficients for the momentum equations
                a[0] = (max(massFlowField_.faceE(i, j, k), 0.) + muField_.faceE(i, j, k)*mesh_.dE(i, j, k)
                        - min(massFlowField_.faceW(i, j, k), 0.) + muField_.faceW(i, j, k)*mesh_.dW(i, j, k)
                        + max(massFlowField_.faceN(i, j, k), 0.) + muField_.faceN(i, j, k)*mesh_.dN(i, j, k)
                        - min(massFlowField_.faceS(i, j, k), 0.) + muField_.faceS(i, j, k)*mesh_.dS(i, j, k)
                        + max(massFlowField_.faceT(i, j, k), 0.) + muField_.faceT(i, j, k)*mesh_.dT(i, j, k)
                        - min(massFlowField_.faceB(i, j, k), 0.) + muField_.faceB(i, j, k)*mesh_.dB(i, j, k)
                        + a0P)/omegaMomentum_;

                a[1] =  min(massFlowField_.faceE(i, j, k), 0.) - muField_.faceE(i, j, k)*mesh_.dE(i, j, k);
                a[2] = -max(massFlowField_.faceW(i, j, k), 0.) - muField_.faceW(i, j, k)*mesh_.dW(i, j, k);
                a[3] =  min(massFlowField_.faceN(i, j, k), 0.) - muField_.faceN(i, j, k)*mesh_.dN(i, j, k);
                a[4] = -max(massFlowField_.faceS(i, j, k), 0.) - muField_.faceS(i, j, k)*mesh_.dS(i, j, k);
                a[5] =  min(massFlowField_.faceT(i, j, k), 0.) - muField_.faceT(i, j, k)*mesh_.dT(i, j, k);
                a[6] = -max(massFlowField_.faceB(i, j, k), 0.) - muField_.faceB(i, j, k)*mesh_.dB(i, j, k);

                bu = a0P*uField0_(i, j, k);
                bu += (1. - omegaMomentum_)*a[0]*uFieldStar_(i, j, k);

                //- Advection/diffusion bcs for the momentum equations
                flowBcs_.uFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, bu);

                for(componentNo = 0; componentNo < 3; ++componentNo)
                {
                    //- Pressure terms for the momentum equations
                    a[7] = mesh_.gE(i, j, k)*mesh_.fAreaNormE(i, j, k)(componentNo)
                            + mesh_.gW(i, j, k)*mesh_.fAreaNormW(i, j, k)(componentNo)
                            + mesh_.gN(i, j, k)*mesh_.fAreaNormN(i, j, k)(componentNo)
                            + mesh_.gS(i, j, k)*mesh_.fAreaNormS(i, j, k)(componentNo)
                            + mesh_.gT(i, j, k)*mesh_.fAreaNormT(i, j, k)(componentNo)
                            + mesh_.gB(i, j, k)*mesh_.fAreaNormB(i, j, k)(componentNo);

                    a[8] = (1. - mesh_.gE(i, j, k))*mesh_.fAreaNormE(i, j, k)(componentNo);
                    a[9] = (1. - mesh_.gW(i, j, k))*mesh_.fAreaNormW(i, j, k)(componentNo);
                    a[10] = (1. - mesh_.gN(i, j, k))*mesh_.fAreaNormN(i, j, k)(componentNo);
                    a[11] = (1. - mesh_.gS(i, j, k))*mesh_.fAreaNormS(i, j, k)(componentNo);
                    a[12] = (1. - mesh_.gT(i, j, k))*mesh_.fAreaNormT(i, j, k)(componentNo);
                    a[13] = (1. - mesh_.gB(i, j, k))*mesh_.fAreaNormB(i, j, k)(componentNo);

                    rowNo = indexMap_(i, j, k, componentNo);
                    cols[0] = rowNo;
                    cols[1] = indexMap_(i + 1, j, k, componentNo);
                    cols[2] = indexMap_(i - 1, j, k, componentNo);
                    cols[3] = indexMap_(i, j + 1, k, componentNo);
                    cols[4] = indexMap_(i, j - 1, k, componentNo);
                    cols[5] = indexMap_(i, j, k + 1, componentNo);
                    cols[6] = indexMap_(i, j, k - 1, componentNo);
                    cols[7] = indexMap_(i, j, k, 3); // Redundant, but not changed for clarity
                    cols[8] = indexMap_(i + 1, j, k, 3);
                    cols[9] = indexMap_(i - 1, j, k, 3);
                    cols[10] = indexMap_(i, j + 1, k, 3);
                    cols[11] = indexMap_(i, j - 1, k, 3);
                    cols[12] = indexMap_(i, j, k + 1, 3);
                    cols[13] = indexMap_(i, j, k - 1, 3);

                    //- Pressure bcs for the momentum equations
                    flowBcs_.pFieldBcs.setImplicitBoundaryCoefficients(i, j, k, &a[7], bu(componentNo));
                    //- Add the momentum terms to the linear system of equations
                    A_[0].setRow(rowNo, 14, cols, a);
                    b_[0].setValue(rowNo, bu(componentNo));
                } // end for componentNo

                dField_(i, j, k) = mesh_.cellVol(i, j, k)/a[0];
            }
        }
    }

    flowBcs_.dFieldBcs.setBoundaries();
    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, dField_);

    //- Continuity equation
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                rowNo = indexMap_(i, j, k, 3);
                bp = 0.;

                for(componentNo = 0; componentNo < 3; ++componentNo)
                {
                    //- Velocity components of the continuity equation
                    a[0] = rhoField_.faceE(i, j, k)*mesh_.gE(i, j, k)*mesh_.fAreaNormE(i, j, k)(componentNo)
                            + rhoField_.faceW(i, j, k)*mesh_.gW(i, j, k)*mesh_.fAreaNormW(i, j, k)(componentNo)
                            + rhoField_.faceN(i, j, k)*mesh_.gN(i, j, k)*mesh_.fAreaNormN(i, j, k)(componentNo)
                            + rhoField_.faceS(i, j, k)*mesh_.gS(i, j, k)*mesh_.fAreaNormS(i, j, k)(componentNo)
                            + rhoField_.faceT(i, j, k)*mesh_.gT(i, j, k)*mesh_.fAreaNormT(i, j, k)(componentNo)
                            + rhoField_.faceB(i, j, k)*mesh_.gB(i, j, k)*mesh_.fAreaNormB(i, j, k)(componentNo);

                    a[1] = rhoField_.faceE(i, j, k)*(1. - mesh_.gE(i, j, k))*mesh_.fAreaNormE(i, j, k)(componentNo);
                    a[2] = rhoField_.faceW(i, j, k)*(1. - mesh_.gW(i, j, k))*mesh_.fAreaNormW(i, j, k)(componentNo);
                    a[3] = rhoField_.faceN(i, j, k)*(1. - mesh_.gN(i, j, k))*mesh_.fAreaNormN(i, j, k)(componentNo);
                    a[4] = rhoField_.faceS(i, j, k)*(1. - mesh_.gS(i, j, k))*mesh_.fAreaNormS(i, j, k)(componentNo);
                    a[5] = rhoField_.faceT(i, j, k)*(1. - mesh_.gT(i, j, k))*mesh_.fAreaNormT(i, j, k)(componentNo);
                    a[6] = rhoField_.faceB(i, j, k)*(1. - mesh_.gB(i, j, k))*mesh_.fAreaNormB(i, j, k)(componentNo);

                    cols[0] = indexMap_(i, j, k, componentNo);
                    cols[1] = indexMap_(i + 1, j, k, componentNo);
                    cols[2] = indexMap_(i - 1, j, k, componentNo);
                    cols[3] = indexMap_(i, j + 1, k, componentNo);
                    cols[4] = indexMap_(i, j - 1, k, componentNo);
                    cols[5] = indexMap_(i, j, k + 1, componentNo);
                    cols[6] = indexMap_(i, j, k - 1, componentNo);

                    bu = Vector3D(0., 0., 0.);
                    flowBcs_.uFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, bu);
                    bp += bu(componentNo);
                    A_[0].setRow(rowNo, 7, cols, a);
                }

                //- Pressure components of the continuity equation, using the Rhie-Chow interpolation
                a[0] = rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*mesh_.dE(i, j, k)
                        + rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*mesh_.dW(i, j, k)
                        + rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*mesh_.dN(i, j, k)
                        + rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*mesh_.dS(i, j, k)
                        + rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*mesh_.dT(i, j, k)
                        + rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*mesh_.dB(i, j, k);

                a[1] = -rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*mesh_.dE(i, j, k);
                a[2] = -rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*mesh_.dW(i, j, k);
                a[3] = -rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*mesh_.dN(i, j, k);
                a[4] = -rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*mesh_.dS(i, j, k);
                a[5] = -rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*mesh_.dT(i, j, k);
                a[6] = -rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*mesh_.dB(i, j, k);

                cols[0] = rowNo;
                cols[1] = indexMap_(i + 1, j, k, 3);
                cols[2] = indexMap_(i - 1, j, k, 3);
                cols[3] = indexMap_(i, j + 1, k, 3);
                cols[4] = indexMap_(i, j - 1, k, 3);
                cols[5] = indexMap_(i, j, k + 1, 3);
                cols[6] = indexMap_(i, j, k - 1, 3);

                bp -= rhoField_.faceE(i, j, k)*dField_.faceE(i, j, k)*dot(gradPField_.faceE(i, j, k), mesh_.fAreaNormE(i, j, k))
                        + rhoField_.faceW(i, j, k)*dField_.faceW(i, j, k)*dot(gradPField_.faceW(i, j, k), mesh_.fAreaNormW(i, j, k))
                        + rhoField_.faceN(i, j, k)*dField_.faceN(i, j, k)*dot(gradPField_.faceN(i, j, k), mesh_.fAreaNormN(i, j, k))
                        + rhoField_.faceS(i, j, k)*dField_.faceS(i, j, k)*dot(gradPField_.faceS(i, j, k), mesh_.fAreaNormS(i, j, k))
                        + rhoField_.faceT(i, j, k)*dField_.faceT(i, j, k)*dot(gradPField_.faceT(i, j, k), mesh_.fAreaNormT(i, j, k))
                        + rhoField_.faceB(i, j, k)*dField_.faceB(i, j, k)*dot(gradPField_.faceB(i, j, k), mesh_.fAreaNormB(i, j, k));

                flowBcs_.pFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, bp);
                A_[0].setRow(rowNo, 7, cols, a);
                b_[0].setValue(rowNo, bp);
            }
        }
    }
    time_.toc();
    Output::print("Coupled", "momentum/continuity matrix assembly time : " + time_.elapsedTime());

    //- Solve the linear system
    time_.tic();
    biCGStabIters_ = A_[0].solve(b_[0], x_[0]);
    time_.toc();
    Output::print("Coupled", "momentum/continuity matrix solution time : " + time_.elapsedTime());

    //- Map the solution back to the mesh
    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                for(componentNo = 0; componentNo < 3; ++componentNo)
                {
                    uField_(i, j, k)(componentNo) = x_[0](indexMap_(i, j, k, componentNo));
                }
                pField_(i, j, k) = x_[0](indexMap_(i, j, k, 3));
            }
        }
    }

    flowBcs_.pFieldBcs.setBoundaries();
    FvScalarScheme::computeCellCenteredGradients(FvScheme::DIVERGENCE_THEOREM, pField_, gradPField_);
    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, gradPField_);

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                hField_(i, j, k) = uField_(i, j, k) + dField_(i, j, k)*gradPField_(i, j, k);
            }
        }
    }

    flowBcs_.hFieldBcs.setBoundaries();
    FvScheme::interpolateInteriorFaces(FvScheme::VOLUME_WEIGHTED, hField_);

    flowBcs_.uFieldBcs.setBoundaries();
    rhieChowInterpolateFaces();
}

void Coupled::rhieChowInterpolateFaces()
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
                    massFlowField_.faceW(i, j, k) = -dot(hField_.faceW(i, j, k), mesh_.fAreaNormW(i, j, k)) - dField_.faceW(i, j, k)*(pField_(i - 1, j, k) - pField_(i, j, k))*mesh_.dW(i, j, k);
                }
                if(j == 0)
                {
                    massFlowField_.faceS(i, j, k) = -dot(hField_.faceS(i, j, k), mesh_.fAreaNormS(i, j, k)) - dField_.faceS(i, j, k)*(pField_(i, j - 1, k) - pField_(i, j, k))*mesh_.dS(i, j, k);
                }
                if(k == 0)
                {
                    massFlowField_.faceB(i, j, k) = -dot(hField_.faceB(i, j, k), mesh_.fAreaNormB(i, j, k)) - dField_.faceB(i, j, k)*(pField_(i, j, k - 1) - pField_(i, j, k))*mesh_.dB(i, j, k);
                }

                massFlowField_.faceE(i, j, k) = dot(hField_.faceE(i, j, k), mesh_.fAreaNormE(i, j, k)) - dField_.faceE(i, j, k)*(pField_(i + 1, j, k) - pField_(i, j, k))*mesh_.dE(i, j, k);
                massFlowField_.faceN(i, j, k) = dot(hField_.faceN(i, j, k), mesh_.fAreaNormN(i, j, k)) - dField_.faceN(i, j, k)*(pField_(i, j + 1, k) - pField_(i, j, k))*mesh_.dN(i, j, k);
                massFlowField_.faceT(i, j, k) = dot(hField_.faceT(i, j, k), mesh_.fAreaNormT(i, j, k)) - dField_.faceT(i, j, k)*(pField_(i, j, k + 1) - pField_(i, j, k))*mesh_.dT(i, j, k);
            }
        }
    }
}

