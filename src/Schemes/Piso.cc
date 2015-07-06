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

void Piso::displayUpdateMessage()
{
    int i, j, k;
    double maxContinuityError = 0.;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) != ACTIVE)
                    continue;

                maxContinuityError = 0.;
            }
        }
    }

    Output::print("Piso", "Momentum prediction total BiCGStab iterations : " + std::to_string(momentumBiCGStabIters_));
    Output::print("Piso", "Pressure correction total BiCGStab iterations : " + std::to_string(pCorrBiCGStabIters_) + "\n");
    Output::print("Piso", "Momentum L1 residual norm     : " + std::to_string(momentumL1Norm_));
    Output::print("Piso", "Momentum L2 residual norm     : " + std::to_string(momentumL2Norm_));
    Output::print("Piso", "Momentum inf residual norm    : " + std::to_string(momentumInfNorm_));
    Output::print("Piso", "Maximum cell continuity error : " + std::to_string(maxContinuityError) + "\n");

    momentumBiCGStabIters_= 0;
    pCorrBiCGStabIters_ = 0;
}

