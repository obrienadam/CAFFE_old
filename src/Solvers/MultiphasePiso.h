#ifndef MULTIPHASE_PISO_H
#define MULTIPHASE_PISO_H

#include "Piso.h"

class MultiphasePiso : public Piso
{
public:

    MultiphasePiso(const Input &input, const HexaFvmMesh &mesh);

    virtual double solve(double timeStep);
    virtual void displayUpdateMessage();

protected:

    void computeAlpha(double timeStep);
    void computeCurvature();

    Field<double> alphaField_;
    Field<Vector3D> gradAlphaField_;
    Field<double> kField_;

    PrimitiveBoundaryCondition<double> alphaFieldBcs_;
};

#endif
