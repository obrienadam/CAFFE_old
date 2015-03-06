#include "Solver.h"

Solver::Solver()
    :
      simTime_(0.),
      itrs_(0.)
{

}

void Solver::initialize(Input &input)
{

}

void Solver::initialize(int nSolutionVariables)
{
    timeDerivatives_.resize(nSolutionVariables);
}

void Solver::reset()
{
    simTime_ = 0.;
    itrs_ = 0;
}

double Solver::simTime()
{
    return simTime_;
}

int Solver::nItrs()
{
    return itrs_;
}
