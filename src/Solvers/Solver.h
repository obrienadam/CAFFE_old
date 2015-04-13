#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

#include "Input.h"
#include "FvScheme.h"

class Solver
{
protected:

    std::vector<double> timeDerivatives_;

    double computeResidualNorm();

public:

    Solver();
    virtual void initialize(Input& input);
    virtual void initialize(int nSolutionVariables);

    virtual double solve(double timeStep, FvScheme& scheme) = 0;
};

#endif
