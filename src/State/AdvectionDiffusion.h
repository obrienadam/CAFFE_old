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

};

#endif // ADVECTIONDIFFUSION_H
