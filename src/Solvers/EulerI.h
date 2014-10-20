#include "Euler.h"

template <class DOMAIN_TYPE, class STATE_TYPE>
Euler<DOMAIN_TYPE, STATE_TYPE>::Euler()
    :
      Solver("Explicit Euler")
{

}

template <class DOMAIN_TYPE, class STATE_TYPE>
void Euler<DOMAIN_TYPE, STATE_TYPE>::initialize(Input &input)
{

    std::ostringstream message;

    Solver::initialize(input);

    Solver::initializeNumOfSteps(1, Solver::domain_.size());

    message << "Initialized solver: " << Solver::solverName_ << "\n"
            << "Fixed time step: " << Solver::timeStep_;

    Output::printLine();

    Output::printToScreen(message.str());

}

template <class DOMAIN_TYPE, class STATE_TYPE>
void Euler<DOMAIN_TYPE, STATE_TYPE>::solve()
{

    int i(0);

    Solver::domain_.computeTimeDerivatives(Solver::timeDerivatives_[0].data());

    for(Solver::itr_ = Solver::begin_; Solver::itr_ != Solver::end_; ++Solver::itr_)
    {

        *(Solver::itr_) += Solver::timeDerivatives_[0][i]*Solver::timeStep_;

        ++i;

    }

}
