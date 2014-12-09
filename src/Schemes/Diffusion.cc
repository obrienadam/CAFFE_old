#include "Diffusion.h"
#include "Output.h"

Diffusion::Diffusion()
{

}

void Diffusion::setMeshPointer(HexaFvmMesh *mesh)
{

    FvScheme::setMeshPointer(mesh);

    phi_ = mesh_->findScalarField("phi");
    mu_ = mesh_->findVectorField("mu");

}
