#include <math.h>

#include "Solver.h"
#include "Output.h"

Solver::Solver()
{

}

double Solver::computeResidualNorm()
{
    double sum = 0.;

    for(int i = 0; i < timeDerivatives_.size(); ++i)
    {
        sum += timeDerivatives_[i]*timeDerivatives_[i];
    }

    if(isnan(sum))
        Output::raiseException("Solver", "computeResidualNorm", "A NaN value has been detected.");

    return sqrt(sum);
}

void Solver::initialize(Input &input)
{

}

void Solver::initialize(int nSolutionVariables)
{
    timeDerivatives_.resize(nSolutionVariables);
}
