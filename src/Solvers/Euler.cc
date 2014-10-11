#include "Euler.h"
#include "Output.h"

Euler::Euler()
    :
      SolverInterface("Explicit Euler")
{



}

void Euler::initialize(Input &input)
{

    timeStep_ = input.inputDoubles["timeStep"];
    timeDerivatives_.resize(1);

    Output::printToScreen("Initialized " + solverName_ + " time integration scheme.");
    Output::printLine();

}

void Euler::advanceSolution(DomainInterface* domain, SchemeInterface* scheme)
{

    timeDerivatives_[0] = domain->computeTimeDerivative(scheme);

    (*domain) += timeDerivatives_[0]*timeStep_;
}
