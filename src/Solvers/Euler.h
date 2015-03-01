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
    void solve(double startTime, double maxTime, double timeStep, int maxItrs, FvScheme &scheme);
};

#endif
