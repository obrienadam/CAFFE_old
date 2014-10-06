#ifndef EULER_FLUID_FIELD_H
#define EULER_FLUID_FIELD_H

#include "ScalarField.h"
#include "VectorField.h"

class EulerFluidField
{

 private:

  //- conserved fields

  ScalarField rho_;
  VectorField mom_;
  ScalarField e_;

  //- primitive fields

  VectorField u_;
  ScalarField mach_;

 public:

};

#endif
