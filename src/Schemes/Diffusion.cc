#include "Diffusion.h"
#include "Output.h"

Diffusion::Diffusion()
{

}

Diffusion::~Diffusion()
{

    // This prevents the destructor from being called on the original object. A little bit shady...

    phiField = Field<double>(0., 0., 0.);

}

void Diffusion::initialize(HexaFvmMesh &mesh, std::string conservedFieldName)
{

    FvScheme::initialize(mesh, conservedFieldName);

    phiField = mesh.findScalarField(conservedFieldName_);

}

double Diffusion::computeFaceFlux(int i, int j, int k, Face face)
{

    Vector3D phiGradient;
    double flux;

}

double Diffusion::computeTimeDerivative(int i, int j, int k)
{



}

void Diffusion::computeSemiDiscreteForm()
{



}
