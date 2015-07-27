#include "SimpleBoundaryCondition.h"

SimpleBoundaryCondition::SimpleBoundaryCondition(const Input &input,
                                                 Field<Vector3D> &uField,
                                                 Field<double> &pField,
                                                 Field<double> &pCorrField,
                                                 Field<double> &rhoField,
                                                 Field<double> &muField,
                                                 Field<Vector3D> &hField,
                                                 Field<double> &dField)
    :
      uFieldBcs(uField),
      pFieldBcs(pField),
      pCorrFieldBcs(pCorrField),
      rhoFieldBcs(rhoField),
      muFieldBcs(muField),
      hFieldBcs(hField),
      dFieldBcs(dField)
{
    using namespace std;

    int i;
    string locationStr[6] = {"east", "west", "north", "south", "top", "bottom"}, typeStr;
    double refPressure;
    Vector3D refVelocity;

    for(i = 0; i < 6; ++i)
    {
        typeStr = input.caseParameters.get<string>("Boundaries." + locationStr[i] + ".type");
        refPressure = input.caseParameters.get<double>("Boundaries." + locationStr[i] + ".refValue");
        refVelocity = stov(input.caseParameters.get<string>("Boundaries." + locationStr[i] + ".refVector"));

        if(typeStr == "inlet")
            types_[i] = INLET;
        else if(typeStr == "outlet")
            types_[i] = OUTLET;
        else if(typeStr == "wall")
            types_[i] = WALL;
        else if(typeStr == "empty")
            types_[i] = EMPTY;
        else if(typeStr == "zeroGradient")
            types_[i] = ZERO_GRADIENT;
        else
            Output::raiseException("SimpleBoundaryCondition", "SimpleBoundaryCondition", "invalid boundary type \"" + typeStr + "\".");

        if(types_[i] == INLET || types_[i] == WALL)
        {
            uFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::FIXED, refVelocity);
            pFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, refPressure);
            pCorrFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0.);
            hFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::FIXED, refVelocity);
        }
        else if(types_[i] == OUTLET)
        {
            uFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::ZERO_GRADIENT, refVelocity);
            pFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::FIXED, refPressure);
            pCorrFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::FIXED, 0.);
            hFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::ZERO_GRADIENT, refVelocity);
        }
        else if(types_[i] == EMPTY)
        {
            uFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::EMPTY, refVelocity);
            pFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::EMPTY, refPressure);
            pCorrFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::EMPTY, 0.);
            hFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::EMPTY, refVelocity);
        }
        else if(types_[i] == ZERO_GRADIENT)
        {
            uFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::ZERO_GRADIENT, refVelocity);
            pFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, refPressure);
            pCorrFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0.);
            hFieldBcs.changeType(i, PrimitiveBoundaryCondition<Vector3D>::ZERO_GRADIENT, refVelocity);
        }

        dFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0.);
        rhoFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0.);
        muFieldBcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0.);
    }
}

void SimpleBoundaryCondition::setImplicitMomentumBoundaryCoefficients(int i, int j, int k, double a[], Vector3D &b)
{
    uFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, b);
}

void SimpleBoundaryCondition::setImplicitPCorrBoundaryCoefficients(int i, int j, int k, double a[], double &b)
{
    pCorrFieldBcs.setImplicitBoundaryCoefficients(i, j, k, a, b);
}
