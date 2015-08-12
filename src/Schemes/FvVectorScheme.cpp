#include "FvVectorScheme.h"
#include "Matrix.h"

void FvVectorScheme::computeCellCenteredGradients(GradientEvaluationMethod method, const Field<Vector3D> &field, Field<Tensor3D> &gradField)
{
    int i, j, k, l, m;
    Matrix A, b, x;
    const HexaFvmMesh &mesh = field.getMesh();

    switch(method)
    {
    case LEAST_SQUARES:

        A.allocate(6, 3);
        b.allocate(6, 1);
        x.allocate(3, 1);

        for(k = 0; k < mesh.nCellsK(); ++k)
        {
            for(j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh.nCellsI(); ++i)
                {
                    b.reallocate(6, 1);

                    //- Construct the least-squares coefficient matrix for the cell
                    A.addVector3DToRow(mesh.rCellE(i, j, k), 0, 0);
                    A.addVector3DToRow(mesh.rCellW(i, j, k), 1, 0);
                    A.addVector3DToRow(mesh.rCellN(i, j, k), 2, 0);
                    A.addVector3DToRow(mesh.rCellS(i, j, k), 3, 0);
                    A.addVector3DToRow(mesh.rCellT(i, j, k), 4, 0);
                    A.addVector3DToRow(mesh.rCellB(i, j, k), 5, 0);

                    for(l = 0; l < 3; ++l)
                    {
                        b(0, 0) = field(i + 1, j, k)(l) - field(i, j, k)(l);
                        b(1, 0) = field(i - 1, j, k)(l) - field(i, j, k)(l);
                        b(2, 0) = field(i, j + 1, k)(l) - field(i, j, k)(l);
                        b(3, 0) = field(i, j - 1, k)(l) - field(i, j, k)(l);
                        b(4, 0) = field(i, j, k + 1)(l) - field(i, j, k)(l);
                        b(5, 0) = field(i, j, k - 1)(l) - field(i, j, k)(l);

                        x = solveLeastSquares(A, b);

                        for(m = 0; m < 3; ++m)
                        {
                            gradField(i, j, k)(l, m) = x(m, 0);
                        }
                    }

                }
            }
        }

        break;
    case DIVERGENCE_THEOREM:

        for(k = 0; k < mesh.nCellsK(); ++k)
        {
            for(j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh.nCellsI(); ++i)
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
    };
}

void FvVectorScheme::extrapolateInteriorFaces(GradientEvaluationMethod method, Field<Vector3D>& field, Field<Tensor3D>& gradField)
{
    using namespace std;

    int i, j, k;
    Vector3D upperLimit, lowerLimit;
    const HexaFvmMesh &mesh = field.getMesh();

    //- Interpolate the interior faces and then use the guessed face values to apply the divergence theorem.
    interpolateInteriorFaces(VOLUME_WEIGHTED, field);
    computeCellCenteredGradients(method, field, gradField);

    for(k = 0; k < mesh.nCellsK(); ++k)
    {
        for(j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh.nCellsI(); ++i)
            {
                if(i < mesh.uCellI())
                {
                    field.faceE(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceE(i, j, k))
                                                + field(i + 1, j, k) + dot(gradField(i + 1, j, k), mesh.rFaceW(i + 1, j, k)));
                }
                if(j < mesh.uCellJ())
                {
                    field.faceN(i, j, k) = 0.5*(field(i, j, k) + dot(gradField(i, j, k), mesh.rFaceN(i, j, k))
                                                + field(i, j + 1, k) + dot(gradField(i, j + 1, k), mesh.rFaceS(i, j + 1, k)));
                }
                if(k < mesh.uCellK())
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
