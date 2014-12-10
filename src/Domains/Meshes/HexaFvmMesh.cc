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

    // Initialize the I-direction faces (normals alligned with in the I-direction)

    nI = nodes_.sizeI();
    nJ = nodes_.sizeJ() - 1;
    nK = nodes_.sizeK() - 1;

    faceCentersI_.allocate(nI, nJ, nK);
    faceNormalsI_.allocate(nI, nJ, nK);
    distanceVectorsI_.allocate(nI, nJ, nK);;
    faceAreasI_.allocate(nI, nJ, nK);
    cellDistancesI_.allocate(nI, nJ, nK);

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

                tmpVec = cellCenters_(i + 1, j, k) - cellCenters_(i, j, k);

                faceCentersI_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsI_(i, j, k) = crossProduct(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                distanceVectorsI_(i, j, k) = tmpVec.unitVector();
                faceAreasI_(i, j, k) = Geometry::computeQuadrilateralArea(tmpPoints);
                cellDistancesI_(i, j, k) = tmpVec.mag();

            } // end for i
        } // end for j
    } // end for k

    // Initialize the J-direction faces (normals aligned with in the J-direction)

    nI = nodes_.sizeI() - 1;
    nJ = nodes_.sizeJ();
    nK = nodes_.sizeK() - 1;

    faceCentersJ_.allocate(nI, nJ, nK);
    faceNormalsJ_.allocate(nI, nJ, nK);
    distanceVectorsJ_.allocate(nI, nJ, nK);
    faceAreasJ_.allocate(nI, nJ, nK);
    cellDistancesJ_.allocate(nI, nJ, nK);

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

                tmpVec = cellCenters_(i, j + 1, k) - cellCenters_(i, j, k);

                faceCentersJ_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsJ_(i, j, k) = crossProduct(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                distanceVectorsJ_(i, j, k) = tmpVec.unitVector();
                faceAreasJ_(i, j, k) = Geometry::computeQuadrilateralArea(tmpPoints);
                cellDistancesJ_(i, j, k) = tmpVec.mag();

            } // end for i
        } // end for j
    } // end for k

    // Initialize the K-direction faces (normals alligned with in the K-direction)

    nI = nodes_.sizeI() - 1;
    nJ = nodes_.sizeJ() - 1;
    nK = nodes_.sizeK();

    faceCentersK_.allocate(nI, nJ, nK);
    faceNormalsK_.allocate(nI, nJ, nK);
    distanceVectorsK_.allocate(nI, nJ, nK);
    faceAreasK_.allocate(nI, nJ, nK);
    cellDistancesK_.allocate(nI, nJ, nK);

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

                tmpVec = cellCenters_(i, j, k + 1) - cellCenters_(i, j, k);

                faceCentersK_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsK_(i, j, k) = crossProduct(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                distanceVectorsK_(i, j, k) = tmpVec.unitVector();
                faceAreasK_(i, j, k) = Geometry::computeQuadrilateralArea(tmpPoints);
                cellDistancesK_(i, j, k) = tmpVec.mag();

            } // end for i
        } // end for j
    } // end for k

    // All fields must now be reallocated

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    for(i = 0; i < scalarFields.size(); ++i)
           scalarFields[i].allocate(nI, nJ, nK);

    for(i = 0; i < vectorFields.size(); ++i)
        vectorFields[i].allocate(nI, nJ, nK);

    nI = faceCentersI_.sizeI();
    nJ = faceCentersI_.sizeJ();
    nK = faceCentersI_.sizeK();

    for(i = 0; i < scalarFluxFieldsI.size(); ++i)
        scalarFluxFieldsI[i].allocate(nI, nJ, nK);

    for(i = 0; i < vectorFluxFieldsI.size(); ++i)
        vectorFluxFieldsI[i].allocate(nI, nJ, nK);

    nI = faceCentersJ_.sizeI();
    nJ = faceCentersJ_.sizeJ();
    nK = faceCentersJ_.sizeK();

    for(i = 0; i < scalarFluxFieldsJ.size(); ++i)
        scalarFluxFieldsJ[i].allocate(nI, nJ, nK);

    for(i = 0; i < vectorFluxFieldsJ.size(); ++i)
        vectorFluxFieldsJ[i].allocate(nI, nJ, nK);

    nI = faceCentersK_.sizeI();
    nJ = faceCentersK_.sizeJ();
    nK = faceCentersK_.sizeK();

    for(i = 0; i < scalarFluxFieldsK.size(); ++i)
        scalarFluxFieldsK[i].allocate(nI, nJ, nK);

    for(i = 0; i < vectorFluxFieldsK.size(); ++i)
        vectorFluxFieldsK[i].allocate(nI, nJ, nK);

    Output::printToScreen("HexaFvmMesh", "Initialization complete.");

}

void HexaFvmMesh::addScalarField(std::string scalarFieldName, int type)
{

    scalarFields.push_back(Field<double>(scalarFieldName, type));
    scalarFields.back().allocate(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK());

    if(type == CONSERVED)
    {

        scalarFluxFieldsI.push_back(Field<double>(scalarFieldName, type));
        scalarFluxFieldsI.back().allocate(faceCentersI_.sizeI(), faceCentersI_.sizeJ(), faceCentersI_.sizeK());

        scalarFluxFieldsJ.push_back(Field<double>(scalarFieldName, type));
        scalarFluxFieldsJ.back().allocate(faceCentersJ_.sizeI(), faceCentersJ_.sizeJ(), faceCentersJ_.sizeK());

        scalarFluxFieldsK.push_back(Field<double>(scalarFieldName, type));
        scalarFluxFieldsK.back().allocate(faceCentersK_.sizeI(), faceCentersK_.sizeJ(), faceCentersK_.sizeK());

    }

    scalarFieldRegistry_[scalarFieldName] = &scalarFields.back();

}

void HexaFvmMesh::addVectorField(std::string vectorFieldName, int type)
{

    vectorFields.push_back(Field<Vector3D>(vectorFieldName));
    vectorFields.back().allocate(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK());

    if(type == CONSERVED)
    {

        vectorFluxFieldsI.push_back(Field<Vector3D>(vectorFieldName, type));
        vectorFluxFieldsI.back().allocate(faceCentersI_.sizeI(), faceCentersI_.sizeJ(), faceCentersI_.sizeK());

        vectorFluxFieldsJ.push_back(Field<Vector3D>(vectorFieldName, type));
        vectorFluxFieldsJ.back().allocate(faceCentersJ_.sizeI(), faceCentersJ_.sizeJ(), faceCentersJ_.sizeK());

        vectorFluxFieldsK.push_back(Field<Vector3D>(vectorFieldName, type));
        vectorFluxFieldsK.back().allocate(faceCentersK_.sizeI(), faceCentersK_.sizeJ(), faceCentersK_.sizeK());

    }

    vectorFieldRegistry_[vectorFieldName] = &vectorFields.back();

}

Field<double>* HexaFvmMesh::findScalarField(std::string fieldName)
{

    if(scalarFieldRegistry_.find(fieldName) == scalarFieldRegistry_.end())
    {

        Output::raiseException("HexaFvmMesh", "findScalarField", "field \"" + fieldName + "\" not found.");

    }

    return scalarFieldRegistry_[fieldName];

}

Field<Vector3D>* HexaFvmMesh::findVectorField(std::string fieldName)
{

    if(vectorFieldRegistry_.find(fieldName) == vectorFieldRegistry_.end())
    {

        Output::raiseException("HexaFvmMesh", "findVectorField", "field \"" + fieldName + "\" not found.");

    }

    return vectorFieldRegistry_[fieldName];

}

Point3D HexaFvmMesh::cellXc(int i, int j, int k)
{

    return cellCenters_(i, j, k);

}

double HexaFvmMesh::cellVol(int i, int j, int k)
{

    return cellVolumes_(i, j, k);

}

Point3D HexaFvmMesh::faceXcE(int i, int j, int k)
{

    return faceCentersI_(i + 1, j, k);

}

Point3D HexaFvmMesh::faceXcW(int i, int j, int k)
{

    return faceCentersI_(i, j, k);

}

Point3D HexaFvmMesh::faceXcN(int i, int j, int k)
{

    return faceCentersJ_(i, j + 1, k);

}

Point3D HexaFvmMesh::faceXcS(int i, int j, int k)
{

    return faceCentersJ_(i, j, k);

}

Point3D HexaFvmMesh::faceXcT(int i, int j, int k)
{


    return faceCentersK_(i, j, k + 1);
}

Point3D HexaFvmMesh::faceXcB(int i, int j, int k)
{

    return faceCentersK_(i, j, k);

}

double HexaFvmMesh::faceAreaE(int i, int j, int k)
{

    return faceAreasI_(i + 1, j, k);

}

double HexaFvmMesh::faceAreaW(int i, int j, int k)
{

    return faceAreasI_(i, j, k);

}

double HexaFvmMesh::faceAreaN(int i, int j, int k)
{

    return faceAreasJ_(i, j + 1, k);

}

double HexaFvmMesh::faceAreaS(int i, int j, int k)
{

    return faceAreasJ_(i, j, k);

}

double HexaFvmMesh::faceAreaT(int i, int j, int k)
{

    return faceAreasK_(i, j, k + 1);

}

double HexaFvmMesh::faceAreaB(int i, int j, int k)
{

    return faceAreasK_(i, j, k);

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
