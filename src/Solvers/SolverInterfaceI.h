#ifndef SOLVER_INTERFACEI_H
#define SOLVER_INTERFACEI_H

#include "SolverInterface.h"

template <class DOMAIN_TYPE, class STATE_TYPE>
void SolverInterface<DOMAIN_TYPE, STATE_TYPE>::initNumOfSteps(int nSteps, int nElements)
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

}

#endif
