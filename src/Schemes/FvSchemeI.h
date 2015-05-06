#include "FvScheme.h"

template <class FIELD>
void FvScheme::interpolateInteriorFaces(FIELD& field, int method)
{
    int i, j, k, uI, uJ, uK;
    HexaFvmMesh& mesh = *meshPtr_;
    double alpha;

    uI = nCellsI_ - 1;
    uJ = nCellsJ_ - 1;
    uK = nCellsK_ - 1;

    switch(method)
    {

    case NON_WEIGHTED:

        alpha = 0.5;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(i < uI)
                    {
                        field.faceE(i, j, k) = alpha*(field(i, j, k) + field(i + 1, j, k));
                    }

                    if(j < uJ)
                    {
                        field.faceN(i, j, k) = alpha*(field(i, j, k) + field(i, j + 1, k));
                    }

                    if(k < uK)
                    {
                        field.faceT(i, j, k) = alpha*(field(i, j, k) + field(i, j, k + 1));
                    }
                }
            }
        }

        break;
    case VOLUME_WEIGHTED:

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(i < uI)
                    {
                        alpha = getAlpha(i, j, k, EAST);
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < uJ)
                    {
                        alpha = getAlpha(i, j, k, NORTH);
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < uK)
                    {
                        alpha = getAlpha(i, j, k, TOP);
                        field.faceT(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j, k + 1);
                    }
                }
            }
        }

        break;
    case DISTANCE_WEIGHTED:

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(i < uI)
                    {
                        alpha = mesh.rFaceMagE(i, j, k)/(mesh.rFaceMagE(i, j, k) + mesh.rFaceMagW(i + 1, j, k));
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < uJ)
                    {
                        alpha = mesh.rFaceMagN(i, j, k)/(mesh.rFaceMagN(i, j, k) + mesh.rFaceMagS(i, j + 1, k));
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < uK)
                    {
                        alpha = mesh.rFaceMagT(i, j, k)/(mesh.rFaceMagT(i, j, k) + mesh.rFaceMagB(i, j, k + 1));
                        field.faceT(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j, k + 1);
                    }
                }
            }
        }

        break;
    default:

        Output::raiseException("FvScheme", "interpolateInteriorFaces", "invalid method selection.");
    };
}
