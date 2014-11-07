#include "SolverInterface.h"

template <class DOMAIN_TYPE, class STATE_TYPE>
void SolverInterface<DOMAIN_TYPE, STATE_TYPE>::initializeNumOfSteps(int nSteps, int nElements)
{

    int k;

    nSteps_ = nSteps;
    nElements_ = nElements;

    timeDerivatives_.resize(nSteps);

    for(k = 0; k < nSteps_; ++k)
    {

        timeDerivatives_[k].resize(nElements_);

    }

}

template <class DOMAIN_TYPE, class STATE_TYPE>
void SolverInterface<DOMAIN_TYPE, STATE_TYPE>::initialize(Input& input)
{

    domain_.initialize(input);

    begin_ = domain_.begin();
    end_ = domain_.end();

    timeStep_ = input.inputDoubles["timeStep"];

}

template <class DOMAIN_TYPE, class STATE_TYPE>
double SolverInterface<DOMAIN_TYPE, STATE_TYPE>::timeStep()
{

    return timeStep_;

}

template <class DOMAIN_TYPE, class STATE_TYPE>
std::string SolverInterface<DOMAIN_TYPE, STATE_TYPE>::timeUnits()
{

    return timeUnits_;

}
