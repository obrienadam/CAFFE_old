#ifndef SIMPLE_BOUNDARY_CONDITION_H
#define SIMPLE_BOUNDARY_CONDITION_H

#include "PrimitiveBoundaryCondition.h"

class SimpleBoundaryCondition
{

public:

    enum Type{INLET, OUTLET, WALL, EMPTY, ZERO_GRADIENT, PARALLEL};

    SimpleBoundaryCondition(const Input &input,
                            Field<Vector3D> &uField,
                            Field<double> &pField,
                            Field<double> &pCorrField,
                            Field<double> &rhoField,
                            Field<double> &muField,
                            Field<Vector3D> &hField,
                            Field<double> &dField);

    Type getTypeEast() const { return types_[0]; }
    Type getTypeWest() const { return types_[1]; }
    Type getTypeNorth() const { return types_[2]; }
    Type getTypeSouth() const { return types_[3]; }
    Type getTypeTop() const { return types_[4]; }
    Type getTypeBottom() const { return types_[5]; }

    void setImplicitMomentumBoundaryCoefficients(int i, int j, int k, double a[], Vector3D &b);
    void setImplicitPCorrBoundaryCoefficients(int i, int j, int k, double a[], double &b);

    PrimitiveBoundaryCondition<Vector3D> uFieldBcs;
    PrimitiveBoundaryCondition<double> pFieldBcs;
    PrimitiveBoundaryCondition<double> pCorrFieldBcs;
    PrimitiveBoundaryCondition<double> rhoFieldBcs;
    PrimitiveBoundaryCondition<double> muFieldBcs;
    PrimitiveBoundaryCondition<Vector3D> hFieldBcs;
    PrimitiveBoundaryCondition<double> dFieldBcs;

private:

    Type types_[6];
};

#endif
