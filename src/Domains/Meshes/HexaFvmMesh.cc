#include "HexaFvmMesh.h"
#include "Geometry.h"

HexaFvmMesh::HexaFvmMesh()
{

}

void HexaFvmMesh::initialize(Input &input)
{

    uint i, j, k;
    Point3D tmpPoints[8];

    // Initialize the mesh nodes

    StructuredMesh::initialize(input);

    // Initialize the cells

    uint nI(nodes_.sizeI() - 1), nJ(nodes_.sizeJ() - 1), nK(nodes_.sizeK() - 1);

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

    // Initialize the J-direction faces (normals alligned with in the J-direction)

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

    // Initialize fields

    for(i = 0; i < scalarFields_.size(); ++i)
    {

        scalarFields_[i].allocate(nI, nJ, nK);

    }

    for(i = 0; i < vectorFields_.size(); ++i)
    {

        vectorFields_[i].allocate(nI, nJ, nK);

    }

}

void HexaFvmMesh::addScalarField(std::string scalarFieldName)
{

    scalarFields_.push_back(Field<double>(scalarFieldName));
    scalarFields_.back().allocate(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK());

}

void HexaFvmMesh::addVectorField(std::string vectorFieldName)
{

    vectorFields_.push_back(Field<Vector3D>(vectorFieldName));
    vectorFields_.back().allocate(cellCenters_.sizeI(), cellCenters_.sizeJ(), cellCenters_.sizeK());

}

void HexaFvmMesh::writeDebug()
{
    using namespace std;

    uint i, j, k, nI, nJ, nK;
    ofstream debugFout;

    debugFout.open((name + "_debug" + ".msh").c_str());
    debugFout << "Cell Positions:\n";

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

    debugFout << "Face Normals I:\n";

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

    debugFout << "Face Normals J:\n";

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

    debugFout << "Face Normals K:\n";

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
