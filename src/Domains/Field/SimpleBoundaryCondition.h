#ifndef SIMPLE_BOUNDARY_CONDITION_H
#define SIMPLE_BOUNDARY_CONDITION_H

#include "FlowBoundaryConditions.h"

class SimpleBoundaryCondition : public FlowBoundaryConditions
{

public:

    SimpleBoundaryCondition(const Input &input,
                            Field<Vector3D> &uField,
                            Field<double> &pField,
                            Field<double> &rhoField,
                            Field<double> &muField,
                            Field<Vector3D> &hField,
                            Field<double> &dField,
                            Field<double> &pCorrField);

    void setImplicitPCorrBoundaryCoefficients(int i, int j, int k, double a[], double &b);

    virtual void setParallelBoundaries(std::shared_ptr< std::array<int, 6> > adjProcNoPtr);

    bool massCorrectionRequiredEast() const;
    bool massCorrectionRequiredWest() const;
    bool massCorrectionRequiredNorth() const;
    bool massCorrectionRequiredSouth() const;
    bool massCorrectionRequiredTop() const;
    bool massCorrectionRequiredBottom() const;

    PrimitiveBoundaryCondition<double> pCorrFieldBcs;
};

#endif
