#include <algorithm>

#include "FvScheme.h"

template <class T>
void FvScheme::interpolateInteriorFaces(InterpolationMethod method, Field<T> &field)
{
    int i, j, k;
    double alpha;
    const HexaFvmMesh &mesh = field.getMesh();

    switch (method)
    {
    case VOLUME_WEIGHTED:

        for(k = 0; k < mesh.nCellsK(); ++k)
        {
            for(j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh.nCellsI(); ++i)
                {
                    if(i < mesh.uCellI())
                    {
                        alpha = mesh.cellVol(i + 1, j, k)/(mesh.cellVol(i + 1, j, k) + mesh.cellVol(i, j, k));

                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < mesh.uCellJ())
                    {
                        alpha = mesh.cellVol(i, j + 1, k)/(mesh.cellVol(i, j + 1, k) + mesh.cellVol(i, j, k));
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < mesh.uCellK())
                    {
                        alpha = mesh.cellVol(i, j, k + 1)/(mesh.cellVol(i, j, k + 1) + mesh.cellVol(i, j, k));
                        field.faceT(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j, k + 1);
                    }
                }
            }
        }

        break;
    case DISTANCE_WEIGHTED:

        for(k = 0; k < mesh.nCellsK(); ++k)
        {
            for(j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh.nCellsI(); ++i)
                {
                    if(i < mesh.uCellI())
                    {
                        alpha = mesh.rFaceMagW(i + 1, j, k)/(mesh.rFaceMagW(i + 1, j, k) + mesh.rFaceMagE(i, j, k));
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < mesh.uCellJ())
                    {
                        alpha = mesh.rFaceMagS(i, j + 1, k)/(mesh.rFaceMagS(i, j + 1, k) + mesh.rFaceMagN(i, j, k));
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < mesh.uCellK())
                    {
                        alpha = mesh.rFaceMagB(i, j, k + 1)/(mesh.rFaceMagB(i, j, k + 1) + mesh.rFaceMagT(i, j, k));
                        field.faceT(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j, k + 1);
                    }
                }
            }
        }

        break;
    case NON_WEIGHTED:

        alpha = 0.5;

        for(k = 0; k < mesh.nCellsK(); ++k)
        {
            for(j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh.nCellsI(); ++i)
                {
                    if(i < mesh.uCellI())
                    {
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < mesh.uCellJ())
                    {
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < mesh.uCellK())
                    {
                        field.faceT(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j, k + 1);
                    }
                }
            }
        }
    };
}
