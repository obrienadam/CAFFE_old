#ifndef ADVECTION_DIFFUSION_FIELD_H
#define ADVECTION_DIFFUSION_FIELD_H

#include "ScalarField.h"
#include "VectorField.h"

class AdvectionDiffusionField
{

 private:

  ScalarField phi_;
  ScalarField alpha_;
  VectorField a_;

 public:

  double& operator()(int i, int j, int k)
  {
    return phi_(i, j, k);
  }

};

#endif
