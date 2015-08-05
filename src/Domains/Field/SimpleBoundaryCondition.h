#ifndef SIMPLE_BOUNDARY_CONDITION_H
#define SIMPLE_BOUNDARY_CONDITION_H

#include "FlowBoundaryConditions.h"

class SimpleBoundaryCondition : public FlowBoundaryConditions
{

public:

    enum Type{INLET, OUTLET, WALL, EMPTY, ZERO_GRADIENT, PARALLEL};

    SimpleBoundaryCondition(const Input &input,
                            Field<Vector3D> &uField,
                            Field<double> &pField,
                            Field<double> &rhoField,
                            Field<double> &muField,
                            Field<Vector3D> &hField,
                            Field<double> &dField,
                            Field<double> &pCorrField);

    void setImplicitPCorrBoundaryCoefficients(int i, int j, int k, double a[], double &b);

    PrimitiveBoundaryCondition<double> pCorrFieldBcs;

private:

    Type types_[6];
};

#endif
