#include "Diffusion.h"
#include "Output.h"

// ************* Private Methods *************

void Diffusion::computeCellCenteredGradients()
{
    int i, j, k, nCellsI(gradPhi_.sizeI()), nCellsJ(gradPhi_.sizeJ()), nCellsK(gradPhi_.sizeK());
    Field<double>& phiField = *phiFieldPtr_;

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                Als_ = cellAls_(i, j, k);

                bls_(0, 0) = phiField(i + 1, j, k) - phiField(i, j, k);
                bls_(1, 0) = phiField(i - 1, j, k) - phiField(i, j, k);
                bls_(2, 0) = phiField(i, j + 1, k) - phiField(i, j, k);
                bls_(3, 0) = phiField(i, j - 1, k) - phiField(i, j, k);
                bls_(4, 0) = phiField(i, j, k + 1) - phiField(i, j, k);
                bls_(5, 0) = phiField(i, j, k - 1) - phiField(i, j, k);

                Als_.solveLeastSquares(bls_);

                gradPhi_(i, j, k).x = bls_(0, 0);
                gradPhi_(i, j, k).y = bls_(1, 0);
                gradPhi_(i, j, k).z = bls_(2, 0);
            } // end for i
        } // end for j
    } // end for k

}

// ************* Public Methods *************

Diffusion::Diffusion()
{

}

Diffusion::~Diffusion()
{

}

void Diffusion::initialize(HexaFvmMesh &mesh, std::string conservedFieldName)
{
    int i, j, k, nCellsI, nCellsJ, nCellsK;

    FvScheme::initialize(mesh, conservedFieldName);

    nCellsI = meshPtr_->nCellsI();
    nCellsJ = meshPtr_->nCellsJ();
    nCellsK = meshPtr_->nCellsK();

    phiFieldPtr_ = &mesh.findScalarField(conservedFieldName_);
    cellAls_.allocate(nCellsI, nCellsJ, nCellsK);

    for(k = 0; k < nCellsK; ++k)
    {
        for(j = 0; j < nCellsJ; ++j)
        {
            for(i = 0; i < nCellsI; ++i)
            {
                // Allocate the least squares matrix

                cellAls_(i, j, k).allocate(6, 3);

                // Populate the least squares matrix

                // East

                if(i != nCellsI - 1)
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceE(i, j, k)*meshPtr_->nesE(i, j, k), 0, 0);
                }
                else
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceE(i, j, k)*meshPtr_->nefE(i, j, k), 0, 0);
                }

                // West

                if(i != 0)
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceW(i, j, k)*meshPtr_->nesW(i, j, k), 1, 0);
                }
                else
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceW(i, j, k)*meshPtr_->nefW(i, j, k), 1, 0);
                }

                // North

                if(j != nCellsJ - 1)
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceN(i, j, k)*meshPtr_->nesN(i, j, k), 2, 0);
                }
                else
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceN(i, j, k)*meshPtr_->nefN(i, j, k), 2, 0);
                }

                // South

                if(j != 0)
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceS(i, j, k)*meshPtr_->nesS(i, j, k), 3, 0);
                }
                else
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceS(i, j, k)*meshPtr_->nefS(i, j, k), 3, 0);
                }

                // Top

                if(k != nCellsK - 1)
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceT(i, j, k)*meshPtr_->nesT(i, j, k), 4, 0);
                }
                else
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceT(i, j, k)*meshPtr_->nefT(i, j, k), 4, 0);
                }

                // Bottom

                if(k != 0)
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceB(i, j, k)*meshPtr_->nesB(i, j, k), 5, 0);
                }
                else
                {
                    cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToFaceDistanceB(i, j, k)*meshPtr_->nefB(i, j, k), 5, 0);
                }

            } // end for i
        } // end for j
    } // end for k

    gradPhi_.allocate(nCellsI, nCellsJ, nCellsK);
    Als_.allocate(6, 3);
    bls_.allocate(6, 1);
}

void Diffusion::discretize()
{

}

void Diffusion::integrate(double timeStep)
{

}
