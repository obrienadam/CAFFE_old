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
    int i, j;

    uField0_ = uField_;

    for(j = 0; j < nInnerIters_; ++j)
    {
        computeMomentum(timeStep);

        for(i = 0; i < nPCorrections_; ++i)
        {
            computePCorr();
            correct();
        }
    }

    return computeContinuityError();
}

void Piso::displayUpdateMessage()
{
    Output::print("Piso", "completed iteration.");
    Output::print("Piso", "BiCGStab iterations: " + std::to_string(biCGStabIters_));
    Output::print("Piso", "Max continuity error: " + std::to_string(continuityError_));
    Output::printLine();
}
