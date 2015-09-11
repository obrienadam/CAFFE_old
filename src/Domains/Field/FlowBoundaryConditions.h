#ifndef FLOW_BOUNDARY_CONDITIONS_H
#define FLOW_BOUNDARY_CONDITIONS_H

#include "PrimitiveBoundaryCondition.h"

class FlowBoundaryConditions
{
public:

    enum Type{INLET, OUTLET, WALL, EMPTY, ZERO_GRADIENT, PARALLEL};

    FlowBoundaryConditions(const Input &input,
                           Field<Vector3D> &uField,
                           Field<double> &pField,
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

    virtual void setParallelBoundaries(std::shared_ptr< std::array<int, 6> > adjProcNoPtr);

    PrimitiveBoundaryCondition<Vector3D> uFieldBcs;
    PrimitiveBoundaryCondition<double> pFieldBcs;
    PrimitiveBoundaryCondition<double> rhoFieldBcs;
    PrimitiveBoundaryCondition<double> muFieldBcs;
    PrimitiveBoundaryCondition<Vector3D> hFieldBcs;
    PrimitiveBoundaryCondition<double> dFieldBcs;

protected:

    Type types_[6];
};

#endif
