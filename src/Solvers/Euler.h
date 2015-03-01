#ifndef EULER_H
#define EULER_H

#include "Solver.h"

class Euler : public Solver
{
private:

public:

    Euler();
    void solve(double maxTime, double timeStep, int maxItrs, FvScheme &scheme);
};

#endif
