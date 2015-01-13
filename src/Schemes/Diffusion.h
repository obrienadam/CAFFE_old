#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "FvScheme.h"
#include "Array3D.h"
#include "Point3D.h"
#include "Vector3D.h"
#include "Matrix.h"

class Diffusion : public FvScheme
{

private:

    //- The matrix containing the least-squares coefficients

    Array3D<Matrix> Als_;
    Matrix xls_;
    Matrix bls_;
    Field<double>* phiFieldPtr_;
    Field<Vector3D> gradPhi_;

    void computeCellCenteredGradients();

public:

    Diffusion();
    ~Diffusion();

    void initialize(HexaFvmMesh &mesh, std::string conservedFieldName);

    double computeFaceFlux(int i, int j, int k, Face face);
    double computeTimeDerivative(int i, int j, int k);
    void computeSemiDiscreteForm();

};

#endif
