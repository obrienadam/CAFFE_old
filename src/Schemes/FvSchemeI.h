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
                    if(cellStatus_(i, j, k) == INACTIVE)
                        continue;

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
                    if(cellStatus_(i, j, k) == INACTIVE)
                        continue;

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
                    if(cellStatus_(i, j, k) == INACTIVE)
                        continue;

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
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                if(i == 0 && field.getWestBoundaryPatch() == ZERO_GRADIENT)
                    field.faceW(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceW(i, j, k));

                if(j == 0 && field.getSouthBoundaryPatch() == ZERO_GRADIENT)
                    field.faceS(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceS(i, j, k));

                if(k == 0 && field.getBottomBoundaryPatch() == ZERO_GRADIENT)
                    field.faceB(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceB(i, j, k));

                if(i == uCellI_ && field.getEastBoundaryPatch() == ZERO_GRADIENT)
                    field.faceE(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceE(i, j, k));
                else if(i < uCellI_)
                {
                    field.faceE(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceE(i, j, k))
                                                + field(i + 1, j, k) + dot(gradField(i + 1, j, k), mesh.rFaceW(i + 1, j, k)));
                }

                if(j == uCellJ_ && field.getNorthBoundaryPatch() == ZERO_GRADIENT)
                    field.faceN(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k));
                else if(j < uCellJ_)
                {
                    field.faceN(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k))
                                                + field(i, j + 1, k) + dot(gradField(i, j + 1, k), mesh.rFaceS(i, j + 1, k)));
                }

                if(k == uCellK_ && field.getTopBoundaryPatch() == ZERO_GRADIENT)
                    field.faceT(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceT(i, j, k));
                else if(k < uCellK_)
                {
                    field.faceT(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceT(i, j, k))
                                                + field(i, j, k + 1) + dot(gradField(i, j, k + 1), mesh.rFaceB(i, j, k + 1)));
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

template <class T, class GRAD_T>
void FvScheme::extrapolateBoundaryFaces(Field<T>& field, Field<GRAD_T>& gradField)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;

    if(field.getEastBoundaryPatch() == ZERO_GRADIENT)
    {
        i = nCellsI_ - 1;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                field.faceE(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceE(i, j, k));
            }
        }
    }

    if(field.getWestBoundaryPatch() == ZERO_GRADIENT)
    {
        i = 0;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(j = 0; j < nCellsJ_; ++j)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                field.faceW(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceW(i, j, k));
            }
        }
    }

    if(field.getNorthBoundaryPatch() == ZERO_GRADIENT)
    {
        j = nCellsJ_ - 1;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                field.faceN(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k));
            }
        }
    }

    if(field.getSouthBoundaryPatch() == ZERO_GRADIENT)
    {
        j = 0;

        for(k = 0; k < nCellsK_; ++k)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                field.faceS(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceS(i, j, k));
            }
        }
    }

    if(field.getTopBoundaryPatch() == ZERO_GRADIENT)
    {
        k = nCellsK_ - 1;

        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                field.faceN(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k));
            }
        }
    }

    if(field.getBottomBoundaryPatch() == ZERO_GRADIENT)
    {
        k = 0;

        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(cellStatus_(i, j, k) == INACTIVE)
                    continue;

                field.faceS(i, j, k) = field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceS(i, j, k));
            }
        }
    }
}

template <class T, class GRAD_T>
void FvScheme::extrapolateAllFaces(Field<T>& field, Field<GRAD_T>& gradField)
{
    extrapolateInteriorFaces(field, gradField);
    extrapolateBoundaryFaces(field, gradField);
}
