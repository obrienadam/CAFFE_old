#include "Euler.h"

Euler::Euler()
{

}

double Euler::solve(double timeStep, FvScheme &scheme)
{
    int i, end;

    if(timeDerivatives_.size() != scheme.nConservedVariables())
        timeDerivatives_.resize(scheme.nConservedVariables());

    end = timeDerivatives_.size();
    scheme.discretize(timeStep, timeDerivatives_);

    for(i = 0; i < end; ++i)
    {
        timeDerivatives_[i] *= timeStep;
    }

    scheme.updateSolution(timeDerivatives_, ADD);

    return computeResidualNorm();
}
