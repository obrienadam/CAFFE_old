#include "LinearAdvection.h"

void LinearAdvection::setMeshPointer(HexaFvmMesh *mesh)
{

    FvScheme::setMeshPointer(mesh);

    phi_ = mesh->findScalarField("phi");
    a_ = mesh->findVectorField("a");

}
