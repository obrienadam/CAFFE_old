#include "FvVectorScheme.h"
#include "Matrix.h"

void FvVectorScheme::computeCellCenteredGradients(GradientEvaluationMethod method, const Field<Vector3D> &field, Field<Tensor3D> &gradField)
{
    const HexaFvmMesh &mesh = field.getMesh();
    Matrix A, b, x;

    switch(method)
    {
    case LEAST_SQUARES:

        A.reallocate(6, 3);
        b.reallocate(6, 1);
        x.reallocate(3, 1);

        for(int k = 0, nK = mesh.nCellsK(); k < nK; ++k)
        {
            for(int j = 0, nJ = mesh.nCellsJ(); j < nJ; ++j)
            {
                for(int i = 0, nI = mesh.nCellsI(); i < nI; ++i)
                {
                    b.reallocate(6, 1);

                    //- Construct the least-squares coefficient matrix for the cell
                    A.addVector3DToRow(mesh.rCellE(i, j, k), 0, 0);
                    A.addVector3DToRow(mesh.rCellW(i, j, k), 1, 0);
                    A.addVector3DToRow(mesh.rCellN(i, j, k), 2, 0);
                    A.addVector3DToRow(mesh.rCellS(i, j, k), 3, 0);
                    A.addVector3DToRow(mesh.rCellT(i, j, k), 4, 0);
                    A.addVector3DToRow(mesh.rCellB(i, j, k), 5, 0);

                    for(int l = 0; l < 3; ++l)
                    {
                        b(0, 0) = field(i + 1, j, k)(l) - field(i, j, k)(l);
                        b(1, 0) = field(i - 1, j, k)(l) - field(i, j, k)(l);
                        b(2, 0) = field(i, j + 1, k)(l) - field(i, j, k)(l);
                        b(3, 0) = field(i, j - 1, k)(l) - field(i, j, k)(l);
                        b(4, 0) = field(i, j, k + 1)(l) - field(i, j, k)(l);
                        b(5, 0) = field(i, j, k - 1)(l) - field(i, j, k)(l);

                        x = solveLeastSquares(A, b);

                        for(int m = 0; m < 3; ++m)
                        {
                            gradField(i, j, k)(l, m) = x(m, 0);
                        }
                    }

                }
            }
        }

        break;

    case DIVERGENCE_THEOREM:

        for(int k = 0, nK = mesh.nCellsK(); k < nK; ++k)
        {
            for(int j = 0, nJ = mesh.nCellsJ(); j < nJ; ++j)
            {
                for(int i = 0, nI = mesh.nCellsI(); i < nI; ++i)
                {
                    gradField(i, j, k) = (tensor(field.faceE(i, j, k), mesh.fAreaNormE(i, j, k))
                                          + tensor(field.faceW(i, j, k), mesh.fAreaNormW(i, j, k))
                                          + tensor(field.faceN(i, j, k), mesh.fAreaNormN(i, j, k))
                                          + tensor(field.faceS(i, j, k), mesh.fAreaNormS(i, j, k))
                                          + tensor(field.faceT(i, j, k), mesh.fAreaNormT(i, j, k))
                                          + tensor(field.faceB(i, j, k), mesh.fAreaNormB(i, j, k)))/mesh.cellVol(i, j, k);
                }
            }
        }

        break;
    }
}

void FvVectorScheme::extrapolateInteriorFaces(GradientEvaluationMethod method, Field<Vector3D>& field, Field<Tensor3D>& gradField)
{
    using namespace std;

    const HexaFvmMesh &mesh = field.getMesh();

    //- Interpolate the interior faces and then use the guessed face values to apply the divergence theorem.
    interpolateInteriorFaces(VOLUME_WEIGHTED, field);
    computeCellCenteredGradients(method, field, gradField);

    for(int k = 0, nK = mesh.nCellsK(); k < nK; ++k)
    {
        for(int j = 0, nJ = mesh.nCellsJ(); j < nJ; ++j)
        {
            for(int i = 0, nI = mesh.nCellsI(); i < nI; ++i)
            {
                bool applyLimitingWest = false, applyLimitingSouth = false, applyLimitingBottom = false;

                if(i == 0 && mesh.westMeshExists())
                {
                    field.faceW(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceW(i, j, k))
                                                + field(i - 1, j, k) + dot(gradField(i - 1, j, k), mesh.rFaceE(i - 1, j, k)));
                    applyLimitingWest = true;
                }

                if(j == 0 && mesh.southMeshExists())
                {
                    field.faceS(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceS(i, j, k))
                                                + field(i, j - 1, k) + dot(gradField(i, j - 1, k), mesh.rFaceN(i, j - 1, k)));
                    applyLimitingSouth = true;
                }

                if(k == 0 && mesh.bottomMeshExists())
                {
                    field.faceB(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceB(i, j, k))
                                                + field(i, j, k - 1) + dot(gradField(i, j, k - 1), mesh.rFaceT(i, j, k - 1)));
                    applyLimitingBottom = true;
                }

                bool applyLimitingEast = false, applyLimitingNorth = false, applyLimitingTop = false;

                if(i < mesh.uCellI() || mesh.eastMeshExists())
                {
                    field.faceE(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceE(i, j, k))
                                                + field(i + 1, j, k) + dot(gradField(i + 1, j, k), mesh.rFaceW(i + 1, j, k)));
                    applyLimitingEast = true;
                }

                if(j < mesh.uCellJ() || mesh.northMeshExists())
                {
                    field.faceN(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k))
                                                + field(i, j + 1, k) + dot(gradField(i, j + 1, k), mesh.rFaceS(i, j + 1, k)));
                    applyLimitingNorth = true;
                }

                if(k < mesh.uCellK() || mesh.topMeshExists())
                {
                    field.faceT(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceT(i, j, k))
                                                + field(i, j, k + 1) + dot(gradField(i, j, k + 1), mesh.rFaceB(i, j, k + 1)));
                    applyLimitingTop = true;
                }

                //- Apply limiting
                Vector3D upperLimit = field.maxNeighbour(i, j, k);
                Vector3D lowerLimit = field.minNeighbour(i, j, k);

                if(applyLimitingEast)
                    field.faceE(i, j, k) = max(min(upperLimit, field.faceE(i, j, k)), lowerLimit);

                if(applyLimitingNorth)
                    field.faceN(i, j, k) = max(min(upperLimit, field.faceN(i, j, k)), lowerLimit);

                if(applyLimitingTop)
                    field.faceT(i, j, k) = max(min(upperLimit, field.faceT(i, j, k)), lowerLimit);

                if(applyLimitingWest)
                    field.faceW(i, j, k) = max(min(upperLimit, field.faceW(i, j, k)), lowerLimit);

                if(applyLimitingSouth)
                    field.faceS(i, j, k) = max(min(upperLimit, field.faceS(i, j, k)), lowerLimit);

                if(applyLimitingBottom)
                    field.faceB(i, j, k) = max(min(upperLimit, field.faceB(i, j, k)), lowerLimit);
            }
        }
    }
}
