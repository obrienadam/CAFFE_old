#ifndef EULER_H
#define EULER_H

#include <vector>

#include "Solver.h"

class Euler : public Solver
{
private:

public:

    Euler();
    double solve(double timeStep, FvScheme &scheme);
};

#endif
