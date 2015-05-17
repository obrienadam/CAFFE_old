#include "PredictorCorrector.h"

PredictorCorrector::PredictorCorrector()
{

}

double PredictorCorrector::solve(double timeStep, FvScheme &scheme)
{
    int i, end;

    if(update_.size() != scheme.nConservedVariables())
    {
        update_.resize(scheme.nConservedVariables());
        original_.resize(scheme.nConservedVariables());
        pTimeDerivatives_.resize(scheme.nConservedVariables());
        cTimeDerivatives_.resize(scheme.nConservedVariables());
        timeDerivatives_.resize(scheme.nConservedVariables());
    }

    scheme.copySolution(original_);

    end = pTimeDerivatives_.size();
    scheme.discretize(timeStep, pTimeDerivatives_);

    for(i = 0; i < end; ++i)
    {
        update_[i] = pTimeDerivatives_[i]*timeStep;
    }

    scheme.updateSolution(update_, ADD);
    scheme.discretize(timeStep, cTimeDerivatives_);

    for(i = 0; i < end; ++i)
    {
        timeDerivatives_[i] = 0.5*(pTimeDerivatives_[i] + cTimeDerivatives_[i]);
        update_[i] = original_[i] + timeStep*timeDerivatives_[i];
    }

    scheme.updateSolution(update_, REPLACE);

    return computeResidualNorm();
}

