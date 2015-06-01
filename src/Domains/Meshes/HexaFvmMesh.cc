/**
 * @file    HexaFvmMesh.cc
 * @author  Adam O'Brien <obrienadam89@gmail.com>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file contains the implementation of methods for class HexaFvmMesh.
 */

#include "HexaFvmMesh.h"
#include "Geometry.h"
#include "Output.h"

// ************* Constructors and Destructors *************

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
    int i, j, k, nI, nJ, nK;
    Vector3D tmpVec;

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    // Allocate the cell to cell distance vectors and distaces

    cellToCellRelativeVectorsI_.allocate(nI - 1, nJ, nK);
    cellToCellUnitVectorsI_.allocate(nI - 1, nJ, nK);
    cellToCellDistancesI_.allocate(nI - 1, nJ, nK);

    cellToCellRelativeVectorsJ_.allocate(nI, nJ - 1, nK);
    cellToCellUnitVectorsJ_.allocate(nI, nJ - 1, nK);
    cellToCellDistancesJ_.allocate(nI, nJ - 1, nK);

    cellToCellRelativeVectorsK_.allocate(nI, nJ, nK - 1);
    cellToCellUnitVectorsK_.allocate(nI, nJ, nK - 1);
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
                    cellToCellRelativeVectorsI_(i, j, k) = tmpVec;
                    cellToCellUnitVectorsI_(i, j, k) = tmpVec.unitVector();
                    cellToCellDistancesI_(i, j, k) = tmpVec.mag();
                }

                if (j < nJ - 1)
                {
                    tmpVec = cellCenters_(i, j + 1, k) - cellCenters_(i, j, k);
                    cellToCellRelativeVectorsJ_(i, j, k) = tmpVec;
                    cellToCellUnitVectorsJ_(i, j, k) = tmpVec.unitVector();
                    cellToCellDistancesJ_(i, j, k) = tmpVec.mag();
                }

                if (k < nK - 1)
                {
                    tmpVec = cellCenters_(i, j, k + 1) - cellCenters_(i, j, k);
                    cellToCellRelativeVectorsK_(i, j, k) = tmpVec;
                    cellToCellUnitVectorsK_(i, j, k) = tmpVec.unitVector();
                    cellToCellDistancesK_(i, j, k) = tmpVec.mag();
                }
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::initializeFaces()
{
    int i, j, k, nI, nJ, nK;
    Point3D tmpPoints[4];
    Vector3D tmpVec;
    double tmpArea;

    // Initialize the I-direction faces (normals alligned with in the I-direction)

    nI = cellCenters_.sizeI() + 1;
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    faceCentersI_.allocate(nI, nJ, nK);
    faceNormalsI_.allocate(nI, nJ, nK);
    faceUnitNormalsI_.allocate(nI, nJ, nK);
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

                tmpVec = cross(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                tmpArea = Geometry::computeQuadrilateralArea(tmpPoints);

                faceCentersI_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsI_(i, j, k) = tmpVec*tmpArea;
                faceUnitNormalsI_(i, j, k) = tmpVec.unitVector();
                faceAreasI_(i, j, k) = tmpArea;
            } // end for i
        } // end for j
    } // end for k

    // Initialize the J-direction faces (normals aligned with in the J-direction)

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ() + 1;
    nK = cellCenters_.sizeK();

    faceCentersJ_.allocate(nI, nJ, nK);
    faceNormalsJ_.allocate(nI, nJ, nK);
    faceUnitNormalsJ_.allocate(nI, nJ, nK);
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

                tmpVec = cross(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                tmpArea = Geometry::computeQuadrilateralArea(tmpPoints);

                faceCentersJ_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsJ_(i, j, k) = tmpVec*tmpArea;
                faceUnitNormalsJ_(i, j, k) = tmpVec.unitVector();
                faceAreasJ_(i, j, k) = tmpArea;
            } // end for i
        } // end for j
    } // end for k

    // Initialize the K-direction faces (normals alligned with in the K-direction)

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK() + 1;

    faceCentersK_.allocate(nI, nJ, nK);
    faceNormalsK_.allocate(nI, nJ, nK);
    faceUnitNormalsK_.allocate(nI, nJ, nK);
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

                tmpVec = cross(tmpPoints[1] - tmpPoints[0], tmpPoints[2] - tmpPoints[0]).unitVector();
                tmpArea = Geometry::computeQuadrilateralArea(tmpPoints);

                faceCentersK_(i, j, k) = Geometry::computeQuadrilateralCentroid(tmpPoints);
                faceNormalsK_(i, j, k) = tmpVec*tmpArea;
                faceUnitNormalsK_(i, j, k) = tmpVec.unitVector();
                faceAreasK_(i, j, k) = tmpArea;
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::initializeCellToFaceParameters()
{
    int i, j, k, nI, nJ, nK;
    Vector3D tmpVec;

    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK();

    cellToFaceRelativeVectorsE_.allocate(nI, nJ, nK);
    cellToFaceUnitVectorsE_.allocate(nI, nJ, nK);
    cellToFaceDistancesE_.allocate(nI, nJ, nK);

    cellToFaceRelativeVectorsW_.allocate(nI, nJ, nK);
    cellToFaceUnitVectorsW_.allocate(nI, nJ, nK);
    cellToFaceDistancesW_.allocate(nI, nJ, nK);

    cellToFaceRelativeVectorsN_.allocate(nI, nJ, nK);
    cellToFaceUnitVectorsN_.allocate(nI, nJ, nK);
    cellToFaceDistancesN_.allocate(nI, nJ, nK);

    cellToFaceRelativeVectorsS_.allocate(nI, nJ, nK);
    cellToFaceUnitVectorsS_.allocate(nI, nJ, nK);
    cellToFaceDistancesS_.allocate(nI, nJ, nK);

    cellToFaceRelativeVectorsT_.allocate(nI, nJ, nK);
    cellToFaceUnitVectorsT_.allocate(nI, nJ, nK);
    cellToFaceDistancesT_.allocate(nI, nJ, nK);

    cellToFaceRelativeVectorsB_.allocate(nI, nJ, nK);
    cellToFaceUnitVectorsB_.allocate(nI, nJ, nK);
    cellToFaceDistancesB_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                // East cell to face parameters

                tmpVec = faceCentersI_(i + 1, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsE_(i, j, k) = tmpVec;
                cellToFaceUnitVectorsE_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesE_(i, j, k) = tmpVec.mag();

                // West cell to face parameters

                tmpVec = faceCentersI_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsW_(i, j, k) = tmpVec;
                cellToFaceUnitVectorsW_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesW_(i, j, k) = tmpVec.mag();

                // North cell to face parameters

                tmpVec = faceCentersJ_(i, j + 1, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsN_(i, j, k) = tmpVec;
                cellToFaceUnitVectorsN_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesN_(i, j, k) = tmpVec.mag();

                // South cell to face parameters

                tmpVec = faceCentersJ_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsS_(i, j, k) = tmpVec;
                cellToFaceUnitVectorsS_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesS_(i, j, k) = tmpVec.mag();

                // Top cell to face parameters

                tmpVec = faceCentersK_(i, j, k + 1) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsT_(i, j, k) = tmpVec;
                cellToFaceUnitVectorsT_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesT_(i, j, k) = tmpVec.mag();

                // Bottom cell to face parameters

                tmpVec = faceCentersK_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsB_(i, j, k) = tmpVec;
                cellToFaceUnitVectorsB_(i, j, k) = tmpVec.unitVector();
                cellToFaceDistancesB_(i, j, k) = tmpVec.mag();
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::initializeGlobalIndexMaps()
{
    // Index maps are used for assembling linear systems. The ordering can be chosen such that the bandwidth is minimized.

    int i, j, k, nI(nodes_.sizeI() - 1), nJ(nodes_.sizeJ() - 1), nK(nodes_.sizeK() - 1);

    rowVectorOrdering_.allocate(nI, nJ, nK);
    columnVectorOrdering_.allocate(nI, nJ, nK);
    layerVectorOrdering_.allocate(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                rowVectorOrdering_(i, j, k) = i + nI*j + nI*nJ*k;
                columnVectorOrdering_(i, j, k) = j + nJ*k + nJ*nK*i;
                layerVectorOrdering_(i, j, k) = k + nK*i + nI*nK*j;
            }// end for i
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
    initializeGlobalIndexMaps();

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

    Output::print("HexaFvmMesh", "Initialization complete.");
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

Point3D HexaFvmMesh::cellXc(int i, int j, int k)
{
    return cellCenters_(i, j, k);
}

Point3D HexaFvmMesh::node(int i, int j, int k, int nodeNo)
{
    switch(nodeNo)
    {
    case BSW:
        return nodes_(i, j, k);
    case BSE:
        return nodes_(i + 1, j, k);
    case BNE:
        return nodes_(i + 1, j + 1, k);
    case BNW:
        return nodes_(i, j + 1, k);
    case TSW:
        return nodes_(i, j, k + 1);
    case TSE:
        return nodes_(i + 1, j, k + 1);
    case TNE:
        return nodes_(i + 1, j + 1, k + 1);
    case TNW:
        return nodes_(i, j + 1, k + 1);
    default:
        Output::raiseException("HexaFvmMesh", "node", "invalid node specified.");
    };

    return Point3D();
}

int HexaFvmMesh::globalIndex(int i, int j, int k, Ordering vectorOrdering)
{
    switch(vectorOrdering)
    {
    case ROW: return rowVectorOrdering_(i, j, k);
    case COLUMN: return columnVectorOrdering_(i, j, k);
    case LAYER: return layerVectorOrdering_(i, j, k);
    };

    return rowVectorOrdering_(i, j, k);
}

Vector3D HexaFvmMesh::rCellE(int i, int j, int k)
{
    if(i == cellToCellRelativeVectorsI_.sizeI())
        return cellToFaceRelativeVectorsE_(i, j, k);

    return cellToCellRelativeVectorsI_(i, j, k);
}

Vector3D HexaFvmMesh::rCellW(int i, int j, int k)
{
    if(i == 0)
        return cellToFaceRelativeVectorsW_(i, j, k);

    return -cellToCellRelativeVectorsI_(i - 1, j, k);
}

Vector3D HexaFvmMesh::rCellN(int i, int j, int k)
{
    if(j == cellToCellRelativeVectorsJ_.sizeJ())
        return cellToFaceRelativeVectorsN_(i, j, k);

    return cellToCellRelativeVectorsJ_(i, j, k);
}

Vector3D HexaFvmMesh::rCellS(int i, int j, int k)
{
    if(j == 0)
        return cellToFaceRelativeVectorsS_(i, j, k);

    return -cellToCellRelativeVectorsJ_(i, j - 1, k);
}

Vector3D HexaFvmMesh::rCellT(int i, int j, int k)
{
    if(k == cellToCellRelativeVectorsK_.sizeK())
        return cellToFaceRelativeVectorsT_(i, j, k);

    return cellToCellRelativeVectorsK_(i, j, k);
}

Vector3D HexaFvmMesh::rCellB(int i, int j, int k)
{
    if(k == 0)
        return cellToFaceRelativeVectorsB_(i, j, k);

    return -cellToCellRelativeVectorsK_(i, j, k - 1);
}

Vector3D HexaFvmMesh::rnCellE(int i, int j, int k)
{
    if(i == cellToCellUnitVectorsI_.sizeI())
        return cellToFaceUnitVectorsE_(i, j, k);

    return cellToCellUnitVectorsI_(i, j, k);
}

Vector3D HexaFvmMesh::rnCellW(int i, int j, int k)
{
    if(i == 0)
        return cellToFaceUnitVectorsW_(i, j, k);

    return -cellToCellUnitVectorsI_(i - 1, j, k);
}

Vector3D HexaFvmMesh::rnCellN(int i, int j, int k)
{
    if(j == cellToCellUnitVectorsJ_.sizeJ())
        return cellToFaceUnitVectorsN_(i, j, k);

    return cellToCellUnitVectorsJ_(i, j, k);
}

Vector3D HexaFvmMesh::rnCellS(int i, int j, int k)
{
    if(j == 0)
        return cellToFaceUnitVectorsS_(i, j, k);

    return -cellToCellUnitVectorsJ_(i, j - 1, k);
}

Vector3D HexaFvmMesh::rnCellT(int i, int j, int k)
{
    if(k == cellToCellUnitVectorsK_.sizeK())
        return cellToFaceUnitVectorsT_(i, j, k);

    return cellToCellUnitVectorsK_(i, j, k);
}

Vector3D HexaFvmMesh::rnCellB(int i, int j, int k)
{
    if(k == 0)
        return cellToFaceUnitVectorsB_(i, j, k);

    return -cellToCellUnitVectorsK_(i, j, k - 1);
}

double HexaFvmMesh::rCellMagE(int i, int j, int k)
{
    if(i == cellToCellDistancesI_.sizeI())
        return cellToFaceDistancesE_(i, j, k);

    return cellToCellDistancesI_(i, j, k);
}

double HexaFvmMesh::rCellMagW(int i, int j, int k)
{
    if(i == 0)
        return cellToFaceDistancesW_(i, j, k);

    return cellToCellDistancesI_(i - 1, j, k);
}

double HexaFvmMesh::rCellMagN(int i, int j, int k)
{
    if(j == cellToCellDistancesJ_.sizeJ())
        return cellToFaceDistancesN_(i, j, k);

    return cellToCellDistancesJ_(i, j, k);
}

double HexaFvmMesh::rCellMagS(int i, int j, int k)
{
    if(j == 0)
        return cellToFaceDistancesS_(i, j, k);

    return cellToCellDistancesJ_(i, j - 1, k);
}

double HexaFvmMesh::rCellMagT(int i, int j, int k)
{
    if(k == cellToCellDistancesK_.sizeK())
        return cellToFaceDistancesT_(i, j, k);

    return cellToCellDistancesK_(i, j, k);
}

double HexaFvmMesh::rCellMagB(int i, int j, int k)
{
    if(k == 0)
        return cellToFaceDistancesB_(i, j, k);

    return cellToCellDistancesK_(i, j, k - 1);
}

void HexaFvmMesh::writeDebug()
{
    using namespace std;

    uint i, j, k;
    ofstream debugFout;

    Output::print("HexaFvmMesh", "writing a debugging file...");

    i = 0;
    j = 0;
    k = 0;

    debugFout.open((name + "_debug" + ".msh").c_str());
    debugFout << "Cell " << i << ", " << j << ", " << k << endl
              << endl
              << "Center: " << cellXc(i, j, k) << endl
              << "Volume: " << cellVol(i, j, k) << endl
              << endl
              << "rCellE: " << rCellE(i, j, k) << endl
              << "rCellW: " << rCellW(i, j, k) << endl
              << "rCellN: " << rCellN(i, j, k) << endl
              << "rCellS: " << rCellS(i, j, k) << endl
              << "rCellT: " << rCellT(i, j, k) << endl
              << "rCellB: " << rCellB(i, j, k) << endl
              << endl
              << "rFaceE: " << rFaceE(i, j, k) << endl
              << "rFaceW: " << rFaceW(i, j, k) << endl
              << "rFaceN: " << rFaceN(i, j, k) << endl
              << "rFaceS: " << rFaceS(i, j, k) << endl
              << "rFaceT: " << rFaceT(i, j, k) << endl
              << "rFaceB: " << rFaceB(i, j, k) << endl
              << endl
              << "fAreaNormE: " << fAreaNormE(i, j, k) << endl
              << "fAreaNormW: " << fAreaNormW(i, j, k) << endl
              << "fAreaNormN: " << fAreaNormN(i, j, k) << endl
              << "fAreaNormS: " << fAreaNormS(i, j, k) << endl
              << "fAreaNormT: " << fAreaNormT(i, j, k) << endl
              << "fAreaNormB: " << fAreaNormB(i, j, k) << endl;

    debugFout.close();

    Output::print("HexaFvmMesh", "finished writing debugging file.");
}

void HexaFvmMesh::writeTec360(double time, std::string directoryName)
{
    using namespace std;

    int i, j, k, varNo, componentNo;
    string component;

    Output::print("HexaFvmMesh", "Writing data to Tec360 ASCII...");

    if(!foutTec360_.is_open())
    {
        foutTec360_.open((directoryName + "/" + name + ".dat").c_str());

        foutTec360_ << "TITLE = \"" << name << "\"" << endl
                    << "VARIABLES = \"x\", \"y\", \"z\", ";

        for(varNo = 0; varNo < scalarFields.size(); ++varNo)
        {
            foutTec360_ << "\"" << scalarFields[varNo].name << "\", ";
        }

        for(varNo = 0; varNo < vectorFields.size(); ++varNo)
        {
            for(componentNo = 0; componentNo < 3; ++componentNo)
            {
                switch(componentNo)
                {
                case 0:
                    component = "x";
                    break;
                case 1:
                    component = "y";
                    break;
                case 2:
                    component = "z";
                    break;
                }

                foutTec360_ << "\"" << vectorFields[varNo].name + "_" + component << "\", ";
            }
        }

        // Special header only for the first zone. This enables the nodes to be shared by subsequent zones
        foutTec360_ << "ZONE T = \"" << name << "Time:" << time << "s\"" << endl
                    << "STRANDID = 1, " << "SOLUTIONTIME = " << time << endl
                    << "I = " << nodes_.sizeI() << ", J = " << nodes_.sizeJ() << ", K = " << nodes_.sizeK() << endl
                    << "DATAPACKING = BLOCK" << endl
                    << "VARLOCATION = ([4-" << 3 + scalarFields.size() + 3*vectorFields.size() << "] = CELLCENTERED)" << endl;

        // Output the mesh data
        for(componentNo = 0; componentNo < 3; ++componentNo)
        {
            for(k = 0; k < nodes_.sizeK(); ++k)
            {
                for(j = 0; j < nodes_.sizeJ(); ++j)
                {
                    for(i = 0; i < nodes_.sizeI(); ++i)
                    {
                        foutTec360_ << nodes_(i, j, k)(componentNo) << " ";
                    }

                    foutTec360_ << endl;
                }
            }
        }
    }
    else
    {
        foutTec360_ << "ZONE T = \"" << name << "Time:" << time << "s\"" << endl
                    << "STRANDID = 1, " << "SOLUTIONTIME = " << time << endl
                    << "I = " << nodes_.sizeI() << ", J = " << nodes_.sizeJ() << ", K = " << nodes_.sizeK() << endl
                    << "DATAPACKING = BLOCK" << endl
                    << "VARSHARELIST = ([1-3] = 1)" << endl
                    << "VARLOCATION = ([4-" << 3 + scalarFields.size() + 3*vectorFields.size() << "] = CELLCENTERED)" << endl;
    }

    // Output the solution data for scalars

    for(varNo = 0; varNo < scalarFields.size(); ++varNo)
    {
        for(k = 0; k < cellCenters_.sizeK(); ++k)
        {
            for(j = 0; j < cellCenters_.sizeJ(); ++j)
            {
                for(i = 0; i < cellCenters_.sizeI(); ++i)
                {
                    foutTec360_ << scalarFields[varNo](i, j, k) << " ";
                }

                foutTec360_ << endl;
            }
        }
    }

    // Output the solution data for vectors

    for(varNo = 0; varNo < vectorFields.size(); ++varNo)
    {
        for(componentNo = 0; componentNo < 3; ++componentNo)
        {
            for(k = 0; k < cellCenters_.sizeK(); ++k)
            {
                for(j = 0; j < cellCenters_.sizeJ(); ++j)
                {
                    for(i = 0; i < cellCenters_.sizeI(); ++i)
                    {
                        foutTec360_ << vectorFields[varNo](i, j, k)(componentNo) << " ";
                    }

                    foutTec360_ << endl;
                }
            }
        }
    }

    Output::print("HexaFvmMesh", "Finished writing data to Tec360 ASCII.");
}
