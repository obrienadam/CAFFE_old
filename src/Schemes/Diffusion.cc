#include "Diffusion.h"
#include "Output.h"

Diffusion::Diffusion()
{

}

void Diffusion::setMeshPointer(HexaFvmMesh *mesh)
{

    int i, j, k;

    FvScheme::setMeshPointer(mesh);

    lsMats_.allocate(mesh->nCellsI(), mesh->nCellsJ(), mesh->nCellsK());

    for(k = 0; k < lsMats_.sizeK(); ++k)
    {

        for(j = 0; j < lsMats_.sizeJ(); ++j)
        {

            for(i = 0; i < lsMats_.sizeI(); ++i)
            {

                lsMats_(i, j, k).allocate(6, 3);

                // Construct the least squares reconstruction coefficient matrices

                // lsMats_(i, j, k).addVector3DToRow(mesh->nesE(i, j, k)*mesh->cellDistanceE(i, j, k), 0, 0);
                // lsMats_(i, j, k).addVector3DToRow(mesh->nesW(i, j, k)*mesh->cellDistanceW(i, j, k), 1, 0);
                // lsMats_(i, j, k).addVector3DToRow(mesh->nesN(i, j, k)*mesh->cellDistanceN(i, j, k), 2, 0);
                // lsMats_(i, j, k).addVector3DToRow(mesh->nesS(i, j, k)*mesh->cellDistanceS(i, j, k), 3, 0);
                // lsMats_(i, j, k).addVector3DToRow(mesh->nesT(i, j, k)*mesh->cellDistanceT(i, j, k), 4, 0);
                // lsMats_(i, j, k).addVector3DToRow(mesh->nesB(i, j, k)*mesh->cellDistanceB(i, j, k), 5, 0);

            } // end for i
        } // end for j
    } // end for k

}
