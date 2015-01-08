#ifndef LINEAR_ADVECTION_H
#define LINEAR_ADVECTION_H

#include <string>

#include "FvScheme.h"
#include "Field.h"
#include "Vector3D.h"

class LinearAdvection : public FvScheme
{

private:

    std::string velocityFieldName_;

public:

    LinearAdvection();

    void initialize(HexaFvmMesh &mesh, std::string conservedFieldName = "phi", std::string velocityFieldName = "a");
    void computeSemiDiscreteForm();

};

#endif
