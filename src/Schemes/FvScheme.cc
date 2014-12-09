#include "FvScheme.h"
#include "Output.h"

FvScheme::FvScheme()
{

}

void FvScheme::setMeshPointer(HexaFvmMesh *mesh)
{

    mesh_ = mesh;

}
