#include "Euler.h"

Euler::Euler()
{

}

void Euler::solve(double timeStep, FvScheme &scheme)
{
    int i, end = timeDerivatives_.size();

    scheme.discretize(timeDerivatives_);

    for(i = 0; i < end; ++i)
    {
        timeDerivatives_[i] *= timeStep;
    }

    scheme.updateSolution(timeDerivatives_);

    simTime_ += timeStep;
    ++itrs_;
}
