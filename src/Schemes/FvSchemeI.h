#include <algorithm>
#include "FvScheme.h"

template <class T>
void FvScheme::interpolateInteriorFaces(Field<T>& field, int method)
{
    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;
    double alpha;

    switch (method)
    {

    case VOLUME_WEIGHTED:

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(i < uCellI_)
                    {
                        alpha = mesh.cellVol(i + 1, j, k)/(mesh.cellVol(i + 1, j, k) + mesh.cellVol(i, j, k));
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < uCellJ_)
                    {
                        alpha = mesh.cellVol(i, j + 1, k)/(mesh.cellVol(i, j + 1, k) + mesh.cellVol(i, j, k));
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < uCellK_)
                    {
                        alpha = mesh.cellVol(i, j, k + 1)/(mesh.cellVol(i, j, k + 1) + mesh.cellVol(i, j, k));
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
                    if(i < uCellI_)
                    {
                        alpha = mesh.rFaceMagW(i + 1, j, k)/(mesh.rFaceMagW(i + 1, j, k) + mesh.rFaceMagE(i, j, k));
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < uCellJ_)
                    {
                        alpha = mesh.rFaceMagS(i, j + 1, k)/(mesh.rFaceMagS(i, j + 1, k) + mesh.rFaceMagN(i, j, k));
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < uCellK_)
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

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                for(i = 0; i < nCellsI_; ++i)
                {
                    if(i < uCellI_)
                    {
                        field.faceE(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i + 1, j, k);
                    }

                    if(j < uCellJ_)
                    {
                        field.faceN(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j + 1, k);
                    }

                    if(k < uCellK_)
                    {
                        field.faceT(i, j, k) = alpha*field(i, j, k) + (1. - alpha)*field(i, j, k + 1);
                    }
                }
            }
        }
    };
}

template <class T, class GRAD_T>
void FvScheme::extrapolateInteriorFaces(Field<T>& field, Field<GRAD_T>& gradField)
{
    using namespace std;

    int i, j, k;
    HexaFvmMesh& mesh = *meshPtr_;
    T upperLimit, lowerLimit;

    //- Interpolate the interior faces and then use the guessed face values to apply the divergence theorem.
    interpolateInteriorFaces(field, VOLUME_WEIGHTED);
    computeCellCenteredGradients(field, gradField, DIVERGENCE_THEOREM);

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {

                //- Extrapolate the faces
                if(i < uCellI_)
                {
                    field.faceE(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceE(i, j, k))
                                                + field(i + 1, j, k) + dot(gradField(i + 1, j, k), mesh.rFaceW(i + 1, j, k)));
                }

                if(j < uCellJ_)
                {
                    field.faceN(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k))
                                                + field(i, j, k) + dot(gradField(i, j + 1, k), mesh.rFaceS(i, j + 1, k)));
                }

                if(k < uCellK_)
                {
                    field.faceT(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceT(i, j, k))
                                                + field(i, j, k) + dot(gradField(i, j, k + 1), mesh.rFaceB(i, j, k + 1)));
                }

                //- Apply limiting
                upperLimit = field.maxNeighbour(i, j, k);
                lowerLimit = field.minNeighbour(i, j, k);

                field.faceE(i, j, k) = max(min(upperLimit, field.faceE(i, j, k)), lowerLimit);
                field.faceN(i, j, k) = max(min(upperLimit, field.faceN(i, j, k)), lowerLimit);
                field.faceT(i, j, k) = max(min(upperLimit, field.faceT(i, j, k)), lowerLimit);
            }
        }
    }
}
