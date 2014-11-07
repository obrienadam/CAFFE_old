#ifndef ADVECTION_DIFFUSION_H
#define ADVECTION_DIFFUSION_H

#include "Vector3D.h"

class AdvectionDiffusion
{

public:

    AdvectionDiffusion(double phi = 298.,
                       double alpha = 1.,
                       Vector3D v = Vector3D(0., 0., 0.));

    double phi, alpha, source;
    Vector3D v;

    AdvectionDiffusion& operator+=(const AdvectionDiffusion& rhs);
    AdvectionDiffusion& operator-=(const AdvectionDiffusion& rhs);
    AdvectionDiffusion& operator*=(double rhs);
    AdvectionDiffusion& operator/=(double rhs);

};

//- Functions

AdvectionDiffusion operator+(AdvectionDiffusion lhs, const AdvectionDiffusion& rhs);
AdvectionDiffusion operator-(AdvectionDiffusion lhs, const AdvectionDiffusion& rhs);
AdvectionDiffusion operator*(AdvectionDiffusion lhs, double rhs);
AdvectionDiffusion operator*(double lhs, AdvectionDiffusion rhs);
AdvectionDiffusion operator/(AdvectionDiffusion lhs, double rhs);

#endif
