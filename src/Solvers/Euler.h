#ifndef EULER_H
#define EULER_H

#include <vector>

#include "Solver.h"

class Euler : public Solver
{
private:

    std::vector<double> timeDerivatives_;

public:

    Euler();
    void solve(double timeStep, FvScheme &scheme);
};

#endif
