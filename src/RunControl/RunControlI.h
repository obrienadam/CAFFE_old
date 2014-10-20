#include "RunControl.h"

template <class DOMAIN_TYPE, class STATE_TYPE>
void RunControl::solverInitialize(SolverInterface<DOMAIN_TYPE, STATE_TYPE>*& solver)
{

    if (input_.inputStrings["solver"] == "Euler")
    {

        solver = new Euler<DOMAIN_TYPE, STATE_TYPE>;

    }
    else
    {

        throw "Unrecognized solver type in RunControl::solverInitialize.";

    }

    solver->initialize(input_);

}
