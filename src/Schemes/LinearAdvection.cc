#include "LinearAdvection.h"

LinearAdvection::LinearAdvection()
{

}

void LinearAdvection::initialize(HexaFvmMesh &mesh, std::string conservedFieldName, std::string velocityFieldName)
{

    FvScheme::initialize(mesh, conservedFieldName);

    velocityFieldName_ = velocityFieldName;

}


void LinearAdvection::computeSemiDiscreteForm()
{

}
