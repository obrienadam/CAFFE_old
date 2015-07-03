#include <iomanip>

#include "Piso.h"

void Piso::discretize(double timeStep, std::vector<double> &timeDerivatives)
{
    Field<Vector3D>& uField = *uFieldPtr_;
    Field<double>& pField = *pFieldPtr_;
    Field<double>& rhoField = *rhoFieldPtr_;
    Field<double>& muField = *muFieldPtr_;
    Field<double>& massFlowField = *massFlowFieldPtr_;
    int i, j;

    storeUField(uField, uField0_);

    for(i = 0; i < maxInnerIters_; ++i)
    {   
        computeMomentum(rhoField, muField, massFlowField, NULL, timeStep, uField, pField);

        for(j = 0; j < 2; ++j)
        {
            computePCorr(rhoField, massFlowField, uField, pField);
            correctContinuity(rhoField, massFlowField, uField, pField);
        }

        std::cout << "\rDTS iteration completion  |      " << (i + 1) << "/" << maxInnerIters_ << std::fixed << std::setprecision(2) << " (" << 100.*(i + 1)/maxInnerIters_ << "%)";
    }

    computeResidual(uField);
}

