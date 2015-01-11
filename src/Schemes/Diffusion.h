#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "FvScheme.h"
#include "Array3D.h"
#include "Point3D.h"
#include "Vector3D.h"

class Diffusion : public FvScheme
{

private:

public:

    Diffusion();
    ~Diffusion();

    Field<double> phiField;

    void initialize(HexaFvmMesh &mesh, std::string conservedFieldName);

    double computeFaceFlux(int i, int j, int k, Face face);
    double computeTimeDerivative(int i, int j, int k);
    void computeSemiDiscreteForm();

};

#endif
