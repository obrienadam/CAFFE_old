#include "HexaFvmMesh.h"
#include "Geometry.h"
#include "Output.h"

// ************* Private Methods *************

void HexaFvmMesh::initializeCells()
{
    int nI(nodes_.sizeI() - 1), nJ(nodes_.sizeJ() - 1), nK(nodes_.sizeK() - 1), i, j, k;
    Point3D tmpPoints[8];

    // Allocate cell centers and their volumes

    cellCenters_.allocate(nI, nJ, nK);
    cellVolumes_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                tmpPoints[0] = nodes_(i, j, k);
                tmpPoints[1] = nodes_(i + 1, j, k);
                tmpPoints[2] = nodes_(i + 1, j + 1, k);
                tmpPoints[3] = nodes_(i, j + 1, k);
                tmpPoints[4] = nodes_(i, j, k + 1);
                tmpPoints[5] = nodes_(i + 1, j, k + 1);
                tmpPoints[6] = nodes_(i + 1, j + 1, k + 1);
                tmpPoints[7] = nodes_(i, j + 1, k + 1);

                cellCenters_(i, j, k) = Geometry::computeHexahedronCentroid(tmpPoints);
                cellVolumes_(i, j, k) = Geometry::computeHexahedronVolume(tmpPoints);
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::initializeCellToCellParameters()
{
    int nI(nodes_.sizeI() - 1), nJ(nodes_.sizeJ() - 1), nK(nodes_.sizeK() - 1), i, j, k;
    Vector3D tmpVec;

    // Allocate the cell to cell distance vectors and distaces

    cellToCellDistanceVectorsI_.allocate(nI - 1, nJ, nK);
    cellToCellDistanceVectorsJ_.allocate(nI, nJ - 1, nK);
    cellToCellDistanceVectorsK_.allocate(nI, nJ, nK - 1);
    cellToCellDistancesI_.allocate(nI - 1, nJ, nK);
    cellToCellDistancesJ_.allocate(nI, nJ - 1, nK);
    cellToCellDistancesK_.allocate(nI, nJ, nK - 1);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                if (i < nI - 1)
                {
                    tmpVec = cellCenters_(i + 1, j, k) - cellCenters_(i, j, k);
                    cellToCellDistanceVectorsI_(i, j, k) = tmpVec.unitVector();
                    cellToCellDistancesI_(i, j, k) = tmpVec.mag();
                }

                if (j < nJ - 1)
                {
                    tmpVec = cellCenters_(i, j + 1, k) - cellCenters_(i, j, k);
                    cellToCellDistanceVectorsJ_(i, j, k) = tmpVec.unitVector();
                    cellToCellDistancesJ_(i, j, k) = tmpVec.mag();
                }

                if (k < nK - 1)
                {
                    tmpVec = cellCenters_(i, j, k + 1) - cellCenters_(i, j, k);
                    cellToCellDistanceVectorsK_(i, j, k) = tmpVec.unitVector();
                    cellToCellDistancesK_(i, j, k) = tmpVec.mag();
                }
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::initializeFaces()
{
    int nI, nJ, nK, i, j, k;
    Point3D tmpPoints[4];

    // Initialize the I-direction faces (normals alligned with in the I-direction)

    nI = nodes_.sizeI();
    nJ = nodes_.sizeJ() - 1;
    nK = nodes_.sizeK() - 1;

    faceCentersI_.allocate(nI, nJ, nK);
    faceNormalsI_.allocate(nI, nJ, nK);
    faceAreasI_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                tmpPoints[0] = nodes_(i, j, k);
                tmpPoints[1] = nodes_(i, j + 1, k);
                tmpPoints[2] = nodes_(i, j + 1, k + 1);
                tmpPoints[3] = nodes_(i, j, k + 1);

                faceCentersI_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsI_(i, j, k) = crossProduct(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                faceAreasI_(i, j, k) = Geometry::computeQuadrilateralArea(tmpPoints);
            } // end for i
        } // end for j
    } // end for k

    // Initialize the J-direction faces (normals aligned with in the J-direction)

    nI = nodes_.sizeI() - 1;
    nJ = nodes_.sizeJ();
    nK = nodes_.sizeK() - 1;

    faceCentersJ_.allocate(nI, nJ, nK);
    faceNormalsJ_.allocate(nI, nJ, nK);
    faceAreasJ_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                tmpPoints[0] = nodes_(i, j, k);
                tmpPoints[1] = nodes_(i, j, k + 1);
                tmpPoints[2] = nodes_(i + 1, j, k + 1);
                tmpPoints[3] = nodes_(i + 1, j, k);

                faceCentersJ_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsJ_(i, j, k) = crossProduct(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                faceAreasJ_(i, j, k) = Geometry::computeQuadrilateralArea(tmpPoints);
            } // end for i
        } // end for j
    } // end for k

    // Initialize the K-direction faces (normals alligned with in the K-direction)

    nI = nodes_.sizeI() - 1;
    nJ = nodes_.sizeJ() - 1;
    nK = nodes_.sizeK();

    faceCentersK_.allocate(nI, nJ, nK);
    faceNormalsK_.allocate(nI, nJ, nK);
    faceAreasK_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                tmpPoints[0] = nodes_(i, j, k);
                tmpPoints[1] = nodes_(i + 1, j, k);
                tmpPoints[2] = nodes_(i + 1, j + 1, k);
                tmpPoints[3] = nodes_(i, j + 1, k);

                faceCentersK_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsK_(i, j, k) = crossProduct(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                faceAreasK_(i, j, k) = Geometry::computeQuadrilateralArea(tmpPoints);
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::initializeCellToFaceParameters()
{
    int nI(nodes_.sizeI() - 1), nJ(nodes_.sizeJ() - 1), nK(nodes_.sizeK() - 1), i, j, k;
    Vector3D tmpVec;

    cellToFaceDistanceVectorsE_.allocate(nI, nJ, nK);
    cellToFaceDistanceVectorsW_.allocate(nI, nJ, nK);
    cellToFaceDistanceVectorsN_.allocate(nI, nJ, nK);
    cellToFaceDistanceVectorsS_.allocate(nI, nJ, nK);
    cellToFaceDistanceVectorsT_.allocate(nI, nJ, nK);
    cellToFaceDistanceVectorsB_.allocate(nI, nJ, nK);

    cellToFaceDistancesE_.allocate(nI, nJ, nK);
    cellToFaceDistancesW_.allocate(nI, nJ, nK);
    cellToFaceDistancesN_.allocate(nI, nJ, nK);
    cellToFaceDistancesS_.allocate(nI, nJ, nK);
    cellToFaceDistancesT_.allocate(nI, nJ, nK);
    cellToFaceDistancesB_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                // East cell to face parameters

                tmpVec = faceCentersI_(i + 1, j, k) - cellCenters_(i, j, k);
                cellToFaceDistanceVectorsE_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesE_(i, j, k) = tmpVec.mag();

                // West cell to face parameters

                tmpVec = faceCentersI_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceDistanceVectorsW_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesW_(i, j, k) = tmpVec.mag();

                // North cell to face parameters

                tmpVec = faceCentersJ_(i, j + 1, k) - cellCenters_(i, j, k);
                cellToFaceDistanceVectorsN_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesN_(i, j, k) = tmpVec.mag();

                // South cell to face parameters

                tmpVec = faceCentersJ_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceDistanceVectorsS_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesS_(i, j, k) = tmpVec.mag();

                // Top cell to face parameters

                tmpVec = faceCentersK_(i, j, k + 1) - cellCenters_(i, j, k);
                cellToFaceDistanceVectorsT_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesT_(i, j, k) = tmpVec.mag();

                // Bottom cell to face parameters

                tmpVec = faceCentersK_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceDistanceVectorsB_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesB_(i, j, k) = tmpVec.mag();
            } // end for i
        } // end for j
    } // end for k
}

// ************* Public Methods *************

void HexaFvmMesh::initialize(Input &input)
{
    uint nI, nJ, nK, i;

    // Initialize the mesh nodes

    StructuredMesh::initialize(input);

    // Initialize the finite volume mesh

    initializeCells();
    initializeCellToCellParameters();
    initializeFaces();
    initializeCellToFaceParameters();

    // All fields must now be reallocated

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    for(i = 0; i < scalarFields.size(); ++i)
    {
        scalarFields[i].allocate(nI, nJ, nK);
    }

    for(i = 0; i < vectorFields.size(); ++i)
    {
        vectorFields[i].allocate(nI, nJ, nK);
    }

    Output::printToScreen("HexaFvmMesh", "Initialization complete.");
}

void HexaFvmMesh::addScalarField(std::string scalarFieldName, int type)
{
    Field<double> newScalarField(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK(), scalarFieldName, type);
    scalarFields.push_back(newScalarField);
}

void HexaFvmMesh::addVectorField(std::string vectorFieldName, int type)
{
    Field<Vector3D> newVectorField(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK(), vectorFieldName, type);
    vectorFields.push_back(newVectorField);
}

Field<double>& HexaFvmMesh::findScalarField(const std::string& fieldName)
{
    int i, end(scalarFields.size());

    for(i = 0; i < end; ++i)
    {
        if(fieldName == scalarFields[i].name)
            return scalarFields[i];
    }

    Output::raiseException("HexaFvmMesh", "findScalarField", "cannot find field \"" + fieldName + "\".");

    // return just to suppress compiler warning

    return scalarFields[end];
}

Field<Vector3D> &HexaFvmMesh::findVectorField(const std::string& fieldName)
{
    int i, end(vectorFields.size());

    for(i = 0; i < end; ++i)
    {
        if(fieldName == vectorFields[i].name)
            return vectorFields[i];
    }

    Output::raiseException("HexaFvmMesh", "findVectorField", "cannot find field \"" + fieldName + "\".");

    // return just to suppress compiler warning

    return vectorFields[end];
}

void HexaFvmMesh::writeDebug()
{
    using namespace std;

    uint i, j, k, nI, nJ, nK;
    ofstream debugFout;

    Output::printToScreen("HexaFvmMesh", "writing a debugging file...");

    debugFout.open((name + "_debug" + ".msh").c_str());
    debugFout << "HexaFvm Mesh Data:\n"
              << "\nCell Positions:\n";

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << cellCenters_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j
    } // end for k

    debugFout << "\nCell Volumes:\n";

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << cellVolumes_(i, j, k) << " ";
            } // end for i

            debugFout << endl;

        } // end for j
    } // end for k

    debugFout << "\nFace Centers I:\n";

    nI = faceCentersI_.sizeI();
    nJ = faceCentersI_.sizeJ();
    nK = faceCentersI_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << faceCentersI_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j
    } // end for k

    debugFout << "\nFace Centers J:\n";

    nI = faceCentersJ_.sizeI();
    nJ = faceCentersJ_.sizeJ();
    nK = faceCentersJ_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << faceCentersJ_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j
    } // end for k

    debugFout << "\nFace Centers K:\n";

    nI = faceCentersK_.sizeI();
    nJ = faceCentersK_.sizeJ();
    nK = faceCentersK_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << faceCentersK_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j
    } // end for k

    debugFout << "\nFace Normals I:\n";

    nI = faceNormalsI_.sizeI();
    nJ = faceNormalsI_.sizeJ();
    nK = faceNormalsI_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << faceNormalsI_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j
    } // end for k

    debugFout << "\nFace Normals J:\n";

    nI = faceNormalsJ_.sizeI();
    nJ = faceNormalsJ_.sizeJ();
    nK = faceNormalsJ_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << faceNormalsJ_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j
    } // end for k

    debugFout << "\nFace Normals K:\n";

    nI = faceNormalsK_.sizeI();
    nJ = faceNormalsK_.sizeJ();
    nK = faceNormalsK_.sizeK();

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                debugFout << faceNormalsK_(i, j, k) << " ";
            } // end for i

            debugFout << endl;
        } // end for j

        debugFout << endl;
    } // end for k

    debugFout.close();

    Output::printToScreen("HexaFvmMesh", "finished writing debugging file.");
}
