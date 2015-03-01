#include "Euler.h"

Euler::Euler()
{

}

void Euler::solve(double maxTime, double timeStep, int maxItrs, FvScheme &scheme)
{
    double time;
    int itrs;

    Output::printToScreen("Euler", "beginning time-marching.");

    for(time = 0., itrs = 0; time < maxTime && itrs < maxItrs; time += timeStep, ++itrs)
    {
        scheme.discretize();
        scheme.integrate(timeStep);
    }

    Output::printToScreen("Euler", "time-marching complete.");
}
