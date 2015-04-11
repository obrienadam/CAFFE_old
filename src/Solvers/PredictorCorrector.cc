#include "PredictorCorrector.h"

PredictorCorrector::PredictorCorrector()
{

}

void PredictorCorrector::solve(double timeStep, FvScheme &scheme)
{
    int i, end;

    if(update_.size() != scheme.nConservedVariables())
    {
        update_.resize(scheme.nConservedVariables());
        original_.resize(scheme.nConservedVariables());
        pTimeDerivatives_.resize(scheme.nConservedVariables());
        cTimeDerivatives_.resize(scheme.nConservedVariables());
    }

    scheme.copySolution(original_);

    end = pTimeDerivatives_.size();
    scheme.discretize(pTimeDerivatives_);

    for(i = 0; i < end; ++i)
    {
        update_[i] = pTimeDerivatives_[i]*timeStep;
    }

    scheme.updateSolution(update_, ADD);
    scheme.discretize(cTimeDerivatives_);

    for(i = 0; i < end; ++i)
    {
        update_[i] = original_[i] + 0.5*timeStep*(pTimeDerivatives_[i] + cTimeDerivatives_[i]);
    }

    scheme.updateSolution(update_, REPLACE);
}

