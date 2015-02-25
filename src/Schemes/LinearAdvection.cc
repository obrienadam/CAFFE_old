#include "LinearAdvection.h"

LinearAdvection::LinearAdvection()
{

}

void LinearAdvection::initialize(HexaFvmMesh &mesh, std::string conservedFieldName, std::string velocityFieldName)
{
    FvScheme::initialize(mesh, conservedFieldName);
    velocityFieldName_ = velocityFieldName;
}

double LinearAdvection::computeFaceFlux(int i, int j, int k, Face face)
{

}

double LinearAdvection::computeTimeDerivative(int i, int j, int k)
{

}

void LinearAdvection::computeSemiDiscreteForm()
{

}
