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

    int i, j, k, nCellsI, nCellsJ, nCellsK;

    FvScheme::initialize(mesh, conservedFieldName);

    nCellsI = meshPtr_->nCellsI();
    nCellsJ = meshPtr_->nCellsJ();
    nCellsK = meshPtr_->nCellsK();

    phiField = mesh.findScalarField(conservedFieldName_);
    Als_.allocate(nCellsI, nCellsJ, nCellsK);

    for(k = 0; k < nCellsK; ++k)
    {

        for(j = 0; j < nCellsJ; ++j)
        {

            for(i = 0; i < nCellsI; ++i)
            {

                // Allocate the least squares matrix

                Als_(i, j, k).allocate(6, 3);

                // Populate the least squares matrix

                // East

                if(i != nCellsI - 1)
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceE(i, j, k)*meshPtr_->nesE(i, j, k), 0, 0);

                }
                else
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceE(i, j, k)*meshPtr_->nefE(i, j, k), 0, 0);

                }

                // West

                if(i != 0)
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceW(i, j, k)*meshPtr_->nesW(i, j, k), 1, 0);

                }
                else
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceW(i, j, k)*meshPtr_->nefW(i, j, k), 1, 0);

                }

                // North

                if(j != nCellsJ - 1)
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceN(i, j, k)*meshPtr_->nesN(i, j, k), 2, 0);

                }
                else
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceN(i, j, k)*meshPtr_->nefN(i, j, k), 2, 0);

                }

                // South

                if(j != 0)
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceS(i, j, k)*meshPtr_->nesS(i, j, k), 3, 0);

                }
                else
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceS(i, j, k)*meshPtr_->nefS(i, j, k), 3, 0);

                }

                // Top

                if(k != nCellsK - 1)
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceT(i, j, k)*meshPtr_->nesT(i, j, k), 4, 0);

                }
                else
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceT(i, j, k)*meshPtr_->nefT(i, j, k), 4, 0);

                }

                // Bottom

                if(k != 0)
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceB(i, j, k)*meshPtr_->nesB(i, j, k), 5, 0);

                }
                else
                {

                    Als_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceB(i, j, k)*meshPtr_->nefB(i, j, k), 5, 0);

                }

            } // end for i
        } // end for j
    } // end for k

    Als_(1, 1, 1).print();

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
