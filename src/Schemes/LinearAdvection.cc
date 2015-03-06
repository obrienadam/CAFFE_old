#include "LinearAdvection.h"

LinearAdvection::LinearAdvection()
{

}

void LinearAdvection::initialize(HexaFvmMesh &mesh, std::string conservedFieldName, std::string velocityFieldName)
{
    FvScheme::initialize(mesh, conservedFieldName);
    velocityFieldName_ = velocityFieldName;
}

int LinearAdvection::nConservedVariables()
{

}

void LinearAdvection::discretize(std::vector<double> &timeDerivatives)
{

}

void LinearAdvection::updateSolution(std::vector<double> &timeDerivatives)
{

}
