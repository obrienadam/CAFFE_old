#ifndef PREDICTOR_CORRECTOR_H
#define PREDICTOR_CORRECTOR_H

#include <vector>

#include "Solver.h"

class PredictorCorrector : public Solver
{
private:

    std::vector<double> update_, original_, pTimeDerivatives_, cTimeDerivatives_;

public:

    PredictorCorrector();

    double solve(double timeStep, FvScheme& scheme);
};

#endif
