#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "FvScheme.h"
#include "Field.h"
#include "Array3D.h"
#include "Vector3D.h"
#include "Matrix.h"

class Diffusion : public FvScheme
{

private:



public:

    Diffusion();
    void computeSemiDiscreteForm();

};

#endif
