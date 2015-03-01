#include "Euler.h"

Euler::Euler()
{

}

void Euler::solve(double startTime, double maxTime, double timeStep, int maxItrs, FvScheme &scheme)
{
    double time;
    int itrs, i, n;

    timeDerivatives_.resize(scheme.nConservedVariables());

    n = timeDerivatives_.size();

    Output::printToScreen("Euler", "beginning time-marching.");

    for(time = startTime, itrs = 0; time < maxTime && itrs < maxItrs; time += timeStep, ++itrs)
    {
        scheme.discretize(timeDerivatives_);
        scheme.updateSolution(timeDerivatives_, timeStep);
    }

    Output::printToScreen("Euler", "time-marching complete.");
}
