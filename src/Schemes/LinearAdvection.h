#ifndef LINEAR_ADVECTION_H
#define LINEAR_ADVECTION_H

#include "FvScheme.h"
#include "Field.h"
#include "Vector3D.h"

class LinearAdvection : public FvScheme
{
private:

    Field<double>* phiFieldPtr_;
    std::string velocityFieldName_;

public:

    LinearAdvection();

    void initialize(Input& input, HexaFvmMesh &mesh, std::string conservedFieldName = "phi", std::string velocityFieldName = "a");
    int nConservedVariables();

    void discretize(std::vector<double>& timeDerivatives);
    void copySolution(std::vector<double>& original);
    void updateSolution(std::vector<double>& timeDerivatives, int method);
};

#endif
