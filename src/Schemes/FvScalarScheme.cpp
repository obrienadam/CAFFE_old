#include "FvScalarScheme.h"
#include "Matrix.h"

void FvScalarScheme::computeCellCenteredGradients(GradientEvaluationMethod method, const Field<double> &field, Field<Vector3D> &gradField)
{
    Matrix A, b;
    const HexaFvmMesh& mesh = field.getMesh();

    switch(method)
    {
    case LEAST_SQUARES:

        A.allocate(6, 3);
        b.allocate(6, 1);

        for(int k = 0, nK = mesh.nCellsK(); k < nK; ++k)
        {
            for(int j = 0, nJ = mesh.nCellsJ(); j < nJ; ++j)
            {
                for(int i = 0, nI = mesh.nCellsI(); i < nI; ++i)
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

        for(int k = 0; k < mesh.nCellsK(); ++k)
        {
            for(int j = 0; j < mesh.nCellsJ(); ++j)
            {
                for(int i = 0; i < mesh.nCellsI(); ++i)
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

    const HexaFvmMesh &mesh = field.getMesh();

    //- Interpolate the interior faces and then use the guessed face values to apply the divergence theorem.
    interpolateInteriorFaces(VOLUME_WEIGHTED, field);
    computeCellCenteredGradients(DIVERGENCE_THEOREM, field, gradField);

    for(int k = 0; k < mesh.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh.nCellsI(); ++i)
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
                double upperLimit = field.maxNeighbour(i, j, k);
                double lowerLimit = field.minNeighbour(i, j, k);

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
