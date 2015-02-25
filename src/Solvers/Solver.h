#ifndef SOLVER_H
#define SOLVER_H

#include "Input.h"
#include "Euler.h"

enum SolverType{EULER};

class Solver
{
protected:

    SolverType currentSolver_;

    uint nIters_;
    Euler euler_;

public:

    Solver();
    void initialize(Input& input);

    void solve();
};

#endif
