#include <iomanip>

#include "Piso.h"

Piso::Piso(const Input &input, const HexaFvmMesh &mesh)
    :
      Simple::Simple(input, mesh)
{
    nPCorrections_ = input.caseParameters.get<int>("Solver.numberOfPressureCorrections");
}

double Piso::solve(double timeStep)
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
    courantNumber_ = computeCourantNumber(timeStep);

    return computeContinuityError();
}

void Piso::displayUpdateMessage()
{
    Output::print("Piso", "completed iteration.");
    Output::print("Piso", "BiCGStab iterations  : " + std::to_string(biCGStabIters_));
    Output::print("Piso", "Max continuity error : " + std::to_string(continuityError_));
    Output::print("Piso", "Max Courant number   : " + std::to_string(courantNumber_));
    Output::printLine();
}
