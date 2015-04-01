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

                cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceE(i, j, k)*meshPtr_->nesE(i, j, k), 0, 0);
                cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceW(i, j, k)*meshPtr_->nesW(i, j, k), 1, 0);
                cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceN(i, j, k)*meshPtr_->nesN(i, j, k), 2, 0);
                cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceS(i, j, k)*meshPtr_->nesS(i, j, k), 3, 0);
                cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceT(i, j, k)*meshPtr_->nesT(i, j, k), 4, 0);
                cellAls_(i, j, k).addVector3DToRow(meshPtr_->cellToCellDistanceB(i, j, k)*meshPtr_->nesB(i, j, k), 5, 0);

            } // end for i
        } // end for j
    } // end for k

    gradPhi_.allocate(nCellsI, nCellsJ, nCellsK);
    Als_.allocate(6, 3);
    bls_.allocate(6, 1);
}

int Diffusion::nConservedVariables()
{
    return phiFieldPtr_->size();
}

void Diffusion::discretize(std::vector<double>& timeDerivatives)
{
    Field<double>& phiField = *phiFieldPtr_;
    HexaFvmMesh& mesh = *meshPtr_;

    computeCellCenteredGradients();
}

void Diffusion::updateSolution(std::vector<double> &timeDerivatives)
{

}
