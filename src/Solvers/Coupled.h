#ifndef COUPLED_H
#define COUPLED_H

#include "Solver.h"
#include "Field.h"
#include "FlowBoundaryConditions.h"

class Coupled : public Solver
{
public:

    Coupled(const Input &input, const HexaFvmMesh &mesh);

    virtual double solve(double timeStep);

protected:

    void computeMomentumAndPressure(double timeStep);
    void rhieChowInterpolateFaces();
    double computeContinuityError();

    //- Primary fields
    Field<Vector3D> uField_;
    Field<double> pField_;
    Field<double> rhoField_;
    Field<double> muField_;

    //- Gradient fields
    Field<Vector3D> gradPField_;

    //- Auxillary fields
    Field<double> dField_;
    Field<Vector3D> hField_;
    Field<double> massFlowField_;

    Field<Vector3D> uField0_;
    Field<Vector3D> uFieldStar_;

    FlowBoundaryConditions flowBcs_;

    int nInnerIters_;
    double omegaMomentum_;
};

#endif
