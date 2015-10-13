#include <iomanip>

#include "MultiphasePiso.h"
#include "InitialConditions.h"

MultiphasePiso::MultiphasePiso(const Input &input, const HexaFvmMesh &mesh)
    :
      Piso::Piso(input, mesh),
      alphaField_(mesh, Field<double>::CONSERVED, "alpha"),
      gradAlphaField_(mesh, Field<Vector3D>::CONSERVED, "gradAlpha"),
      kField_(mesh, Field<double>::CONSERVED, "k"),
      alphaFieldBcs_(alphaField_)
{
    mesh_.addArray3DToTecplotOutput(alphaField_.name, alphaField_.cellData());
    mesh_.addArray3DToTecplotOutput(gradAlphaField_.name, gradAlphaField_.cellData());
    mesh_.addArray3DToTecplotOutput(kField_.name, kField_.cellData());

    for(int boundaryNo = 0; boundaryNo < 6; ++boundaryNo)
        alphaFieldBcs_.changeType(boundaryNo, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0.);

    alphaFieldBcs_.setParallelBoundaries(mesh_.getAdjProcNoPtr());

    InitialConditions initialConditions;
    initialConditions.setInitialConditions(alphaField_);
}

double MultiphasePiso::solve(double timeStep)
{
    uField0_ = uField_;

    for(int innerIterNo = 0; innerIterNo < nInnerIters_; ++innerIterNo)
    {
        computeMomentum(timeStep);

        for(int pCorrNo = 0; pCorrNo < nPCorrections_; ++pCorrNo)
        {
            computePCorr();
            correct();
        }
    }
    computeAlpha(timeStep);
    computeCurvature();

    courantNumber_ = computeCourantNumber(timeStep);

    return computeContinuityError();
}

void MultiphasePiso::displayUpdateMessage()
{
    Output::print("MultiphasePiso", "completed iteration.");
    Output::print("MultiphasePiso", "BiCGStab iterations  : " + std::to_string(biCGStabIters_));
    Output::print("MultiphasePiso", "Max continuity error : " + std::to_string(continuityError_));
    Output::print("MultiphasePiso", "Max Courant number   : " + std::to_string(courantNumber_));
    Output::printLine();
}

void MultiphasePiso::computeAlpha(double timeStep)
{
    using namespace std;

    int cols[7];
    double a0P, a[7];

    for(int k = 0, nCellsK = mesh_.nCellsK(); k < nCellsK; ++k)
    {
        for(int j = 0, nCellsJ = mesh_.nCellsJ(); j < nCellsJ; ++j)
        {
            for(int i = 0, nCellsI = mesh_.nCellsI(); i < nCellsI; ++i)
            {
                if(!mesh_.iMap.isActive(i, j, k))
                    continue;
                else if(Solver::solutionType_ == UNSTEADY)
                    a0P = mesh_.cellVol(i, j, k)/timeStep;
                else
                    a0P = 0.;

                a[1] =  min(dot(uField_.faceE(i, j, k), mesh_.fAreaNormE(i, j, k)), 0.);
                a[2] =  min(dot(uField_.faceW(i, j, k), mesh_.fAreaNormW(i, j, k)), 0.);
                a[3] =  min(dot(uField_.faceN(i, j, k), mesh_.fAreaNormN(i, j, k)), 0.);
                a[4] =  min(dot(uField_.faceS(i, j, k), mesh_.fAreaNormS(i, j, k)), 0.);
                a[5] =  min(dot(uField_.faceT(i, j, k), mesh_.fAreaNormT(i, j, k)), 0.);
                a[6] =  min(dot(uField_.faceB(i, j, k), mesh_.fAreaNormB(i, j, k)), 0.);

                a[0] = max(dot(uField_.faceE(i, j, k), mesh_.fAreaNormE(i, j, k)), 0.)
                        + max(dot(uField_.faceW(i, j, k), mesh_.fAreaNormW(i, j, k)), 0.)
                        + max(dot(uField_.faceN(i, j, k), mesh_.fAreaNormN(i, j, k)), 0.)
                        + max(dot(uField_.faceS(i, j, k), mesh_.fAreaNormS(i, j, k)), 0.)
                        + max(dot(uField_.faceT(i, j, k), mesh_.fAreaNormT(i, j, k)), 0.)
                        + max(dot(uField_.faceB(i, j, k), mesh_.fAreaNormB(i, j, k)), 0.)
                        + a0P;

                double b = 0.;

                alphaFieldBcs_.setImplicitBoundaryCoefficients(i, j, k, a, b);

                int rowNo = mesh_.iMap(i, j, k, 0);
                cols[0] = mesh_.iMap(i, j, k, 0);
                cols[1] = mesh_.iMap(i + 1, j, k, 0);
                cols[2] = mesh_.iMap(i - 1, j, k, 0);
                cols[3] = mesh_.iMap(i, j + 1, k, 0);
                cols[4] = mesh_.iMap(i, j - 1, k, 0);
                cols[5] = mesh_.iMap(i, j, k + 1, 0);
                cols[6] = mesh_.iMap(i, j, k - 1, 0);

                A_[0].setRow(rowNo, 7, cols, a);
                b_[0].setValue(rowNo, b);
            }
        }
    }
}

void MultiphasePiso::computeCurvature()
{

}
