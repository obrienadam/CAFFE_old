#ifndef EULER_H
#define EULER_H

#include "SolverInterface.h"

template <class DOMAIN_TYPE, class STATE_TYPE>
class Euler : public SolverInterface<DOMAIN_TYPE, STATE_TYPE>
{

private:

    typedef SolverInterface<DOMAIN_TYPE, STATE_TYPE> Solver;

public:

    Euler();

    void initialize(Input &input);
    void solve();

};

#include "EulerI.h"

#endif
