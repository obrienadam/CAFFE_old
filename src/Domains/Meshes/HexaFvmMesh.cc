#include "HexaFvmMesh.h"
#include "Geometry.h"
#include "Output.h"

void HexaFvmMesh::initialize(Input &input)
{

    uint i, j, k, nI, nJ, nK;
    Point3D tmpPoints[8];
    Vector3D tmpVec;

    // Initialize the mesh nodes

    StructuredMesh::initialize(input);

    // Initialize the cells

    nI = nodes_.sizeI() - 1;
    nJ = nodes_.sizeJ() - 1;
    nK = nodes_.sizeK() - 1;

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

    // All fields must now be reallocated

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    for(i = 0; i < scalarFields.size(); ++i)
    {

        scalarFields[i].resize(nI, nJ, nK);

    }

    for(i = 0; i < vectorFields.size(); ++i)
    {

        vectorFields[i].resize(nI, nJ, nK);

    }

    Output::printToScreen("HexaFvmMesh", "Initialization complete.");

}

void HexaFvmMesh::addScalarField(std::string scalarFieldName, int type)
{

    scalarFields.push_back(Field<double>(scalarFieldName, type));
    scalarFields.back().allocate(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK());
    scalarFieldRegistry_[scalarFieldName] = &scalarFields.back();

}

void HexaFvmMesh::addVectorField(std::string vectorFieldName, int type)
{

    vectorFields.push_back(Field<Vector3D>(vectorFieldName));
    vectorFields.back().allocate(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK());
    vectorFieldRegistry_[vectorFieldName] = &vectorFields.back();

}

Field<double>& HexaFvmMesh::findScalarField(std::string fieldName)
{

    if(scalarFieldRegistry_.find(fieldName) == scalarFieldRegistry_.end())
    {

        Output::raiseException("HexaFvmMesh", "findScalarField", "field \"" + fieldName + "\" not found.");

    }

    return *scalarFieldRegistry_[fieldName];

}

Field<Vector3D>& HexaFvmMesh::findVectorField(std::string fieldName)
{

    if(vectorFieldRegistry_.find(fieldName) == vectorFieldRegistry_.end())
    {

        Output::raiseException("HexaFvmMesh", "findVectorField", "field \"" + fieldName + "\" not found.");

    }

    return *vectorFieldRegistry_[fieldName];

}

void HexaFvmMesh::writeDebug()
{
    using namespace std;

    uint i, j, k, nI, nJ, nK;
    ofstream debugFout;

    debugFout.open((name + "_debug" + ".msh").c_str());

    debugFout << "HexaFvm Mesh Data:\n";

    debugFout << "\nCell Positions:\n";

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
    } // end for k

    debugFout.close();

}
