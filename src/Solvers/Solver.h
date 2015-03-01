#ifndef SOLVER_H
#define SOLVER_H

#include "Input.h"
#include "FvScheme.h"

class Solver
{
protected:

public:

    Solver();
    void initialize(Input& input);

    virtual void solve(double startTime, double maxTime, double timeStep, int maxItrs, FvScheme& scheme) = 0;
};

#endif
