#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

#include "Input.h"
#include "FvScheme.h"

class Solver
{
protected:

    std::vector<double> timeDerivatives_;

    double simTime_;
    int itrs_;

public:

    Solver();
    virtual void initialize(Input& input);
    virtual void initialize(int nSolutionVariables);

    virtual void solve(double timeStep, FvScheme& scheme) = 0;
    void reset();
    double simTime();
    int nItrs();
};

#endif
