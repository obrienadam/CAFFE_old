#include "Solver.h"

Solver::Solver()
    :
      nIters_(0)
{



}

void Solver::initialize(Input &input)
{

    if(input.inputStrings["solver"] == "Euler")
    {

        currentSolver_ = EULER;

    }

}

void Solver::solve()
{



}
