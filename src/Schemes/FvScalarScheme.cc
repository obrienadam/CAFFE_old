#include "FvScalarScheme.h"
#include "Matrix.h"

void FvScalarScheme::computeCellCenteredGradients(GradientEvaluationMethod method, const Field<double> &field, Field<Vector3D> &gradField)
{
    int i, j, k;
    Matrix A, b;
    const HexaFvmMesh& mesh = field.getMesh();

    switch(method)
    {
    case LEAST_SQUARES:

        A.allocate(6, 3);
        b.allocate(6, 1);

        for(k = 0; k < mesh.nCellsK(); ++k)
        {
            for(j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(i = 0; i < mesh.nCellsI(); ++i)
                {
                    b.reallocate(6, 1);

                    A.addVector3DToRow(mesh.rCellE(i, j, k), 0, 0);
                    A.addVector3DToRow(mesh.rCellW(i, j, k), 1, 0);
                    A.addVector3DToRow(mesh.rCellN(i, j, k), 2, 0);
                    A.addVector3DToRow(mesh.rCellS(i, j, k), 3, 0);
                    A.addVector3DToRow(mesh.rCellT(i, j, k), 4, 0);
                    A.addVector3DToRow(mesh.rCellB(i, j, k), 5, 0);

                    b(0, 0) = field(i + 1, j, k) - field(i, j, k);
                    b(1, 0) = field(i - 1, j, k) - field(i, j, k);
                    b(2, 0) = field(i, j + 1, k) - field(i, j, k);
                    b(3, 0) = field(i, j - 1, k) - field(i, j, k);
                    b(4, 0) = field(i, j, k + 1) - field(i, j, k);
                    b(5, 0) = field(i, j, k - 1) - field(i, j, k);

                    A.solveLeastSquares(b);

                    gradField(i, j, k).x = b(0, 0);
                    gradField(i, j, k).y = b(1, 0);
                    gradField(i, j, k).z = b(2, 0);
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
                    gradField(i, j, k) = (field.faceE(i, j, k)*mesh.fAreaNormE(i, j, k)
                                          + field.faceW(i, j, k)*mesh.fAreaNormW(i, j, k)
                                          + field.faceN(i, j, k)*mesh.fAreaNormN(i, j, k)
                                          + field.faceS(i, j, k)*mesh.fAreaNormS(i, j, k)
                                          + field.faceT(i, j, k)*mesh.fAreaNormT(i, j, k)
                                          + field.faceB(i, j, k)*mesh.fAreaNormB(i, j, k))/mesh.cellVol(i, j, k);
                }
            }
        }

        break;
    };
}

void FvScalarScheme::extrapolateInteriorFaces(GradientEvaluationMethod method, Field<double>& field, Field<Vector3D>& gradField)
{
    using namespace std;

    int i, j, k;
    double upperLimit, lowerLimit;
    const HexaFvmMesh &mesh = field.getMesh();

    //- Interpolate the interior faces and then use the guessed face values to apply the divergence theorem.
    interpolateInteriorFaces(VOLUME_WEIGHTED, field);
    computeCellCenteredGradients(DIVERGENCE_THEOREM, field, gradField);

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
