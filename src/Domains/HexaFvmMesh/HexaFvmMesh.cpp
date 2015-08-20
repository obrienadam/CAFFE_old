/**
 * @file    HexaFvmMesh.cpp
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

HexaFvmMesh::HexaFvmMesh()
    :
      eastBoundaryMeshPtr_(NULL),
      westBoundaryMeshPtr_(NULL),
      northBoundaryMeshPtr_(NULL),
      southBoundaryMeshPtr_(NULL),
      topBoundaryMeshPtr_(NULL),
      bottomBoundaryMeshPtr_(NULL)
{

}

HexaFvmMesh::HexaFvmMesh(const HexaFvmMesh &other)
    :
      HexaFvmMesh()
{

}

// ************* Public Methods *************

void HexaFvmMesh::initialize(const std::string &filename)
{
    // Initialize the mesh nodes
    StructuredMesh::initialize(filename);
    Output::print("HexaFvmMesh", "initializing hexahedral finite-volume mesh...");
    initializeCellsAndFaces();
    Output::print("HexaFvmMesh", "initialization of hexahedral finite-volume mesh complete.");
}

void HexaFvmMesh::initialize(const Array3D<Point3D> &nodes)
{
    // Initialize the mesh nodes
    StructuredMesh::initialize(nodes);
    initializeCellsAndFaces();
}

void HexaFvmMesh::initialize(const std::vector<Point3D> &vertices, int nCellsI, int nCellsJ, int nCellsK)
{
    StructuredMesh::initialize(vertices, nCellsI + 1, nCellsJ + 1, nCellsK + 1);
    initializeCellsAndFaces();
}

void HexaFvmMesh::initializeCartesianMesh(double xLength, double yLength, double zLength, int nCellsI, int nCellsJ, int nCellsK)
{
    StructuredMesh::initializeCartesianMesh(xLength, yLength, zLength, nCellsI + 1, nCellsJ + 1, nCellsK + 1);
    initializeCellsAndFaces();
}

void HexaFvmMesh::addBoundaryMesh(const HexaFvmMesh &mesh, Direction relativeLocation)
{
    switch(relativeLocation)
    {
    case EAST:
        eastBoundaryMeshPtr_ = &mesh;
        break;
    case WEST:
        westBoundaryMeshPtr_ = &mesh;
        break;
    case NORTH:
        northBoundaryMeshPtr_ = &mesh;
        break;
    case SOUTH:
        southBoundaryMeshPtr_ = &mesh;
        break;
    case TOP:
        topBoundaryMeshPtr_ = &mesh;
        break;
    case BOTTOM:
        bottomBoundaryMeshPtr_ = &mesh;
        break;
    };

    //- Re-compute mesh metrics to account for the presence of the boundary mesh
    computeMeshMetrics();
}

bool HexaFvmMesh::eastMeshExists() const
{
    if(eastBoundaryMeshPtr_ == NULL)
        return false;
    else
        return true;
}

bool HexaFvmMesh::westMeshExists() const
{
    if(westBoundaryMeshPtr_ == NULL)
        return false;
    else
        return true;
}

bool HexaFvmMesh::northMeshExists() const
{
    if(northBoundaryMeshPtr_ == NULL)
        return false;
    else
        return true;
}

bool HexaFvmMesh::southMeshExists() const
{
    if(southBoundaryMeshPtr_ == NULL)
        return false;
    else
        return true;
}

bool HexaFvmMesh::topMeshExists() const
{
    if(topBoundaryMeshPtr_ == NULL)
        return false;
    else
        return true;
}

bool HexaFvmMesh::bottomMeshExists() const
{
    if(bottomBoundaryMeshPtr_ == NULL)
        return false;
    else
        return true;
}

std::string HexaFvmMesh::meshStats()
{
    std::string structuredMeshStats = StructuredMesh::meshStats();
    std::ostringstream stats;

    stats << structuredMeshStats
          << "Cells in I direction -> " << cellCenters_.sizeI() << "\n"
          << "Cells in J direction -> " << cellCenters_.sizeJ() << "\n"
          << "Cells in K direction -> " << cellCenters_.sizeK() << "\n"
          << "Cells total -> " << cellCenters_.size() << "\n";

    return stats.str();
}

Point3D HexaFvmMesh::cellXc(int i, int j, int k) const
{
    if(i < 0 && westBoundaryMeshPtr_ != NULL)
        return westBoundaryMeshPtr_->cellXc(westBoundaryMeshPtr_->uCellI() + i + 1, j, k);
    else if(i > uCellI() && eastBoundaryMeshPtr_ != NULL)
        return eastBoundaryMeshPtr_->cellXc(i - uCellI() - 1, j, k);

    if(j < 0 && southBoundaryMeshPtr_ != NULL)
        return southBoundaryMeshPtr_->cellXc(i, southBoundaryMeshPtr_->uCellJ() + j + 1, k);
    else if(j > uCellJ() && northBoundaryMeshPtr_ != NULL)
        return northBoundaryMeshPtr_->cellXc(i, j - uCellJ() - 1, k);

    if(k < 0 && bottomBoundaryMeshPtr_ != NULL)
        return bottomBoundaryMeshPtr_->cellXc(i, j, bottomBoundaryMeshPtr_->uCellK() + k + 1);
    else if(k > uCellK() && topBoundaryMeshPtr_ != NULL)
        return topBoundaryMeshPtr_->cellXc(i, j, k - uCellK() - 1);

    return cellCenters_(i, j, k);
}

double HexaFvmMesh::cellVol(int i, int j, int k) const
{
    if(i < 0 && westBoundaryMeshPtr_ != NULL)
        return westBoundaryMeshPtr_->cellVol(westBoundaryMeshPtr_->uCellI() + i + 1, j, k);
    else if(i > uCellI() && eastBoundaryMeshPtr_ != NULL)
        return eastBoundaryMeshPtr_->cellVol(i - uCellI() - 1, j, k);

    if(j < 0 && southBoundaryMeshPtr_ != NULL)
        return southBoundaryMeshPtr_->cellVol(i, southBoundaryMeshPtr_->uCellJ() + j + 1, k);
    else if(j > uCellJ() && northBoundaryMeshPtr_ != NULL)
        return northBoundaryMeshPtr_->cellVol(i, j - uCellJ() - 1, k);

    if(k < 0 && bottomBoundaryMeshPtr_ != NULL)
        return bottomBoundaryMeshPtr_->cellVol(i, j, bottomBoundaryMeshPtr_->uCellK() + k + 1);
    else if(k > uCellK() && topBoundaryMeshPtr_ != NULL)
        return topBoundaryMeshPtr_->cellVol(i, j, k - uCellK() - 1);

    return cellVolumes_(i, j, k);
}

Point3D HexaFvmMesh::node(int i, int j, int k, int nodeNo) const
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

Vector3D HexaFvmMesh::rCellE(int i, int j, int k) const
{
    if(i == uCellI())
    {
        if(eastBoundaryMeshPtr_ != NULL)
            return eastBoundaryMeshPtr_->cellXc(0, j, k) - cellXc(i, j, k);
        else
            return cellToFaceRelativeVectorsE_(i, j, k);
    }

    return cellToCellRelativeVectorsI_(i, j, k);
}

Vector3D HexaFvmMesh::rCellW(int i, int j, int k) const
{
    if(i == 0)
    {
        if(westBoundaryMeshPtr_ != NULL)
            return westBoundaryMeshPtr_->cellXc(westBoundaryMeshPtr_->uCellI(), j, k) - cellXc(i, j, k);
        else
            return cellToFaceRelativeVectorsW_(i, j, k);
    }

    return -cellToCellRelativeVectorsI_(i - 1, j, k);
}

Vector3D HexaFvmMesh::rCellN(int i, int j, int k) const
{
    if(j == uCellJ())
    {
        if(northBoundaryMeshPtr_ != NULL)
            return northBoundaryMeshPtr_->cellXc(i, 0, k) - cellXc(i, j, k);
        else
            return cellToFaceRelativeVectorsN_(i, j, k);
    }

    return cellToCellRelativeVectorsJ_(i, j, k);
}

Vector3D HexaFvmMesh::rCellS(int i, int j, int k) const
{
    if(j == 0)
    {
        if(southBoundaryMeshPtr_ != NULL)
            return southBoundaryMeshPtr_->cellXc(i, southBoundaryMeshPtr_->uCellJ(), k) - cellXc(i, j, k);
        else
            return cellToFaceRelativeVectorsS_(i, j, k);
    }

    return -cellToCellRelativeVectorsJ_(i, j - 1, k);
}

Vector3D HexaFvmMesh::rCellT(int i, int j, int k) const
{
    if(k == uCellK())
    {
        if(topBoundaryMeshPtr_ != NULL)
            return topBoundaryMeshPtr_->cellXc(i, j, 0) - cellXc(i, j, k);
        else
            return cellToFaceRelativeVectorsT_(i, j, k);
    }

    return cellToCellRelativeVectorsK_(i, j, k);
}

Vector3D HexaFvmMesh::rCellB(int i, int j, int k) const
{
    if(k == 0)
    {
        if(bottomBoundaryMeshPtr_ != NULL)
            return bottomBoundaryMeshPtr_->cellXc(i, j, bottomBoundaryMeshPtr_->uCellK()) - cellXc(i, j, k);
        else
            return cellToFaceRelativeVectorsB_(i, j, k);
    }

    return -cellToCellRelativeVectorsK_(i, j, k - 1);
}

void HexaFvmMesh::locateCell(const Point3D &point, int &ii, int &jj, int &kk) const
{
    int i, j, k;
    int nI = cellCenters_.sizeI(), nJ = cellCenters_.sizeJ(), nK = cellCenters_.sizeK();
    Point3D tmpPoints[8];

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

                if(Geometry::isInsideHexahedron(point, tmpPoints))
                {
                    ii = i;
                    jj = j;
                    kk = k;
                    return;
                }
            }
        }
    }

    Output::raiseException("HexaFvmMesh", "locateCell", "the specified point \"" + std::to_string(point) + "\" is not inside the domain.");
}

void HexaFvmMesh::locateEnclosingCells(const Point3D &point, int ii[], int jj[], int kk[]) const
{
    int i, j, k, i0, j0, k0;
    Point3D tmpPoints[8];

    locateCell(point, i0, j0, k0);

    for(k = k0 - 1; k <= k0; ++k)
    {
        for(j = j0 - 1; j <= j0; ++j)
        {
            for(i = i0 - 1; i <= i0; ++i)
            {
                tmpPoints[0] = cellCenters_(i, j, k);
                tmpPoints[1] = cellCenters_(i + 1, j, k);
                tmpPoints[2] = cellCenters_(i + 1, j + 1, k);
                tmpPoints[3] = cellCenters_(i, j + 1, k);
                tmpPoints[4] = cellCenters_(i, j, k + 1);
                tmpPoints[5] = cellCenters_(i + 1, j, k + 1);
                tmpPoints[6] = cellCenters_(i + 1, j + 1, k + 1);
                tmpPoints[7] = cellCenters_(i, j + 1, k + 1);

                if(Geometry::isInsideHexahedron(point, tmpPoints))
                {
                    ii[0] = i; jj[0] = j; kk[0] = k;
                    ii[1] = i + 1; jj[1] = j; kk[1] = k;
                    ii[2] = i + 1; jj[2] = j + 1; kk[2] = k;
                    ii[3] = i; jj[3] = j + 1; kk[3] = k;
                    ii[4] = i; jj[4] = j; kk[4] = k + 1;
                    ii[5] = i + 1; jj[5] = j; kk[5] = k + 1;
                    ii[6] = i + 1; jj[6] = j + 1; kk[6] = k + 1;
                    ii[7] = i; jj[7] = j + 1; kk[7] = k + 1;
                    return;
                }
            }
        }
    }

    std::cout << i0 << " " << j0 << " " << k0 << std::endl;
    Output::raiseException("HexaFvmMesh", "locateEnclosingCells", "a problem occurred when trying to locate the enclosing cells for point \"" + std::to_string(point) + "\".");
}

void HexaFvmMesh::writeDebug() const
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

void HexaFvmMesh::addArray3DToTecplotOutput(std::string name, const Array3D<double> *array3DPtr) const
{
    using namespace std;

    scalarVariablePtrs_.push_back(pair<string, const Array3D<double>*>(name, array3DPtr));
}

void HexaFvmMesh::addArray3DToTecplotOutput(std::string name, const Array3D<Vector3D> *array3DPtr) const
{
    using namespace std;

    vectorVariablePtrs_.push_back(pair<string, const Array3D<Vector3D>*>(name, array3DPtr));
}

void HexaFvmMesh::writeTec360(double time, const std::string &directory)
{
    using namespace std;

    int i, j, k, varNo, componentNo, nVariables = 3 + scalarVariablePtrs_.size() + 3*vectorVariablePtrs_.size();
    string component;

    Output::print("HexaFvmMesh", "Writing data to Tec360 ASCII...");

    if(!foutTec360_.is_open())
    {
        foutTec360_.open(directory + name + ".dat");

        foutTec360_ << "TITLE = \"" << name << "\"" << endl
                    << "VARIABLES = \"x\", \"y\", \"z\", ";

        for(varNo = 0; varNo < scalarVariablePtrs_.size(); ++varNo)
        {
            foutTec360_ << "\"" << scalarVariablePtrs_[varNo].first << "\", ";
        }

        for(varNo = 0; varNo < vectorVariablePtrs_.size(); ++varNo)
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

                foutTec360_ << "\"" << vectorVariablePtrs_[varNo].first + "_" + component << "\", ";
            }
        }

        // Special header only for the first zone. This enables the nodes to be shared by subsequent zones
        foutTec360_ << "ZONE T = \"" << name << "Time:" << time << "s\"" << endl
                    << "STRANDID = 1, " << "SOLUTIONTIME = " << time << endl
                    << "I = " << nodes_.sizeI() << ", J = " << nodes_.sizeJ() << ", K = " << nodes_.sizeK() << endl
                    << "DATAPACKING = BLOCK" << endl;

        if(nVariables == 4)
        {
            foutTec360_ << "VARLOCATION = ([4] = CELLCENTERED)" << endl;
        }
        else if(nVariables > 4)
        {
            foutTec360_ << "VARLOCATION = ([4-" << nVariables << "] = CELLCENTERED)" << endl;
        }

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
                    << "VARSHARELIST = ([1-3] = 1)" << endl;

        if(nVariables == 4)
        {
            foutTec360_ << "VARLOCATION = ([4] = CELLCENTERED)" << endl;
        }
        else if(nVariables > 4)
        {
            foutTec360_ << "VARLOCATION = ([4-" << nVariables << "] = CELLCENTERED)" << endl;
        }
    }

    // Output the solution data for scalars
    for(varNo = 0; varNo < scalarVariablePtrs_.size(); ++varNo)
    {
        for(k = 0; k < cellCenters_.sizeK(); ++k)
        {
            for(j = 0; j < cellCenters_.sizeJ(); ++j)
            {
                for(i = 0; i < cellCenters_.sizeI(); ++i)
                {
                    foutTec360_ << scalarVariablePtrs_[varNo].second->operator()(i, j, k) << " ";
                }

                foutTec360_ << endl;
            }
        }
    }

    // Output the solution data for vectors

    for(varNo = 0; varNo < vectorVariablePtrs_.size(); ++varNo)
    {
        for(componentNo = 0; componentNo < 3; ++componentNo)
        {
            for(k = 0; k < cellCenters_.sizeK(); ++k)
            {
                for(j = 0; j < cellCenters_.sizeJ(); ++j)
                {
                    for(i = 0; i < cellCenters_.sizeI(); ++i)
                    {
                        foutTec360_ << vectorVariablePtrs_[varNo].second->operator()(i, j, k)(componentNo) << " ";
                    }

                    foutTec360_ << endl;
                }
            }
        }
    }

    Output::print("HexaFvmMesh", "Finished writing data to Tec360 ASCII.");
}

// ************* Private Methods *************

void HexaFvmMesh::initializeCellsAndFaces()
{
    // Initialize the finite volume mesh
    initializeCells();
    initializeCellToCellParameters();
    initializeFaces();
    initializeCellToFaceParameters();
    computeMeshMetrics();
    iMap.initialize(nCellsI(), nCellsJ(), nCellsK());
}

void HexaFvmMesh::initializeCells()
{
    int nI(nodes_.sizeI() - 1), nJ(nodes_.sizeJ() - 1), nK(nodes_.sizeK() - 1), i, j, k;
    Point3D tmpPoints[8];

    // Allocate cell centers and their volumes
    cellCenters_.resize(nI, nJ, nK);
    cellVolumes_.resize(nI, nJ, nK);

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
    cellToCellRelativeVectorsI_.resize(nI - 1, nJ, nK);
    cellToCellRelativeVectorsJ_.resize(nI, nJ - 1, nK);
    cellToCellRelativeVectorsK_.resize(nI, nJ, nK - 1);

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
                }

                if (j < nJ - 1)
                {
                    tmpVec = cellCenters_(i, j + 1, k) - cellCenters_(i, j, k);
                    cellToCellRelativeVectorsJ_(i, j, k) = tmpVec;
                }

                if (k < nK - 1)
                {
                    tmpVec = cellCenters_(i, j, k + 1) - cellCenters_(i, j, k);
                    cellToCellRelativeVectorsK_(i, j, k) = tmpVec;
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

    faceCentersI_.resize(nI, nJ, nK);
    faceNormalsI_.resize(nI, nJ, nK);

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
            } // end for i
        } // end for j
    } // end for k

    // Initialize the J-direction faces (normals aligned with in the J-direction)
    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ() + 1;
    nK = cellCenters_.sizeK();

    faceCentersJ_.resize(nI, nJ, nK);
    faceNormalsJ_.resize(nI, nJ, nK);

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
            } // end for i
        } // end for j
    } // end for k

    // Initialize the K-direction faces (normals alligned with in the K-direction)
    nI = cellCenters_.sizeI();
    nJ = cellCenters_.sizeJ();
    nK = cellCenters_.sizeK() + 1;

    faceCentersK_.resize(nI, nJ, nK);
    faceNormalsK_.resize(nI, nJ, nK);

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

    cellToFaceRelativeVectorsE_.resize(nI, nJ, nK);
    cellToFaceRelativeVectorsW_.resize(nI, nJ, nK);
    cellToFaceRelativeVectorsN_.resize(nI, nJ, nK);
    cellToFaceRelativeVectorsS_.resize(nI, nJ, nK);
    cellToFaceRelativeVectorsT_.resize(nI, nJ, nK);
    cellToFaceRelativeVectorsB_.resize(nI, nJ, nK);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                // East cell to face parameters
                tmpVec = faceCentersI_(i + 1, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsE_(i, j, k) = tmpVec;

                // West cell to face parameters
                tmpVec = faceCentersI_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsW_(i, j, k) = tmpVec;

                // North cell to face parameters
                tmpVec = faceCentersJ_(i, j + 1, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsN_(i, j, k) = tmpVec;

                // South cell to face parameters
                tmpVec = faceCentersJ_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsS_(i, j, k) = tmpVec;

                // Top cell to face parameters
                tmpVec = faceCentersK_(i, j, k + 1) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsT_(i, j, k) = tmpVec;

                // Bottom cell to face parameters
                tmpVec = faceCentersK_(i, j, k) - cellCenters_(i, j, k);
                cellToFaceRelativeVectorsB_(i, j, k) = tmpVec;
            } // end for i
        } // end for j
    } // end for k
}

void HexaFvmMesh::computeMeshMetrics()
{
    int i, j, k;

    dE_.resize(nCellsI(), nCellsJ(), nCellsK());
    dW_.resize(nCellsI(), nCellsJ(), nCellsK());
    dN_.resize(nCellsI(), nCellsJ(), nCellsK());
    dS_.resize(nCellsI(), nCellsJ(), nCellsK());
    dT_.resize(nCellsI(), nCellsJ(), nCellsK());
    dB_.resize(nCellsI(), nCellsJ(), nCellsK());

    cE_.resize(nCellsI(), nCellsJ(), nCellsK());
    cW_.resize(nCellsI(), nCellsJ(), nCellsK());
    cN_.resize(nCellsI(), nCellsJ(), nCellsK());
    cS_.resize(nCellsI(), nCellsJ(), nCellsK());
    cT_.resize(nCellsI(), nCellsJ(), nCellsK());
    cB_.resize(nCellsI(), nCellsJ(), nCellsK());

    gE_.resize(nCellsI(), nCellsJ(), nCellsK());
    gW_.resize(nCellsI(), nCellsJ(), nCellsK());
    gN_.resize(nCellsI(), nCellsJ(), nCellsK());
    gS_.resize(nCellsI(), nCellsJ(), nCellsK());
    gT_.resize(nCellsI(), nCellsJ(), nCellsK());
    gB_.resize(nCellsI(), nCellsJ(), nCellsK());

    for(k = 0; k < nCellsK(); ++k)
    {
        for(j = 0; j < nCellsJ(); ++j)
        {
            for(i = 0; i < nCellsI(); ++i)
            {
                dE_(i, j, k) = dot(fAreaNormE(i, j, k), fAreaNormE(i, j, k))/dot(fAreaNormE(i, j, k), rCellE(i, j, k));
                dW_(i, j, k) = dot(fAreaNormW(i, j, k), fAreaNormW(i, j, k))/dot(fAreaNormW(i, j, k), rCellW(i, j, k));
                dN_(i, j, k) = dot(fAreaNormN(i, j, k), fAreaNormN(i, j, k))/dot(fAreaNormN(i, j, k), rCellN(i, j, k));
                dS_(i, j, k) = dot(fAreaNormS(i, j, k), fAreaNormS(i, j, k))/dot(fAreaNormS(i, j, k), rCellS(i, j, k));
                dT_(i, j, k) = dot(fAreaNormT(i, j, k), fAreaNormT(i, j, k))/dot(fAreaNormT(i, j, k), rCellT(i, j, k));
                dB_(i, j, k) = dot(fAreaNormB(i, j, k), fAreaNormB(i, j, k))/dot(fAreaNormB(i, j, k), rCellB(i, j, k));

                cE_(i, j, k) = fAreaNormE(i, j, k) - rCellE(i, j, k)*dot(fAreaNormE(i, j, k), fAreaNormE(i, j, k))/dot(fAreaNormE(i, j, k), rCellE(i, j, k));
                cW_(i, j, k) = fAreaNormW(i, j, k) - rCellW(i, j, k)*dot(fAreaNormW(i, j, k), fAreaNormW(i, j, k))/dot(fAreaNormW(i, j, k), rCellW(i, j, k));
                cN_(i, j, k) = fAreaNormN(i, j, k) - rCellN(i, j, k)*dot(fAreaNormN(i, j, k), fAreaNormN(i, j, k))/dot(fAreaNormN(i, j, k), rCellN(i, j, k));
                cS_(i, j, k) = fAreaNormS(i, j, k) - rCellS(i, j, k)*dot(fAreaNormS(i, j, k), fAreaNormS(i, j, k))/dot(fAreaNormS(i, j, k), rCellS(i, j, k));
                cT_(i, j, k) = fAreaNormT(i, j, k) - rCellT(i, j, k)*dot(fAreaNormT(i, j, k), fAreaNormT(i, j, k))/dot(fAreaNormT(i, j, k), rCellT(i, j, k));
                cB_(i, j, k) = fAreaNormB(i, j, k) - rCellB(i, j, k)*dot(fAreaNormB(i, j, k), fAreaNormB(i, j, k))/dot(fAreaNormB(i, j, k), rCellB(i, j, k));

                if(i < uCellI() || eastBoundaryMeshPtr_ != NULL)
                    gE_(i, j, k) = cellVol(i + 1, j, k)/(cellVol(i + 1, j, k) + cellVol(i, j, k));
                else
                    gE_(i, j, k) = 0.;

                if(i > 0 || westBoundaryMeshPtr_ != NULL)
                    gW_(i, j, k) = cellVol(i - 1, j, k)/(cellVol(i - 1, j, k) + cellVol(i, j, k));
                else
                    gW_(i, j, k) = 0.;

                if(j < uCellJ() || northBoundaryMeshPtr_ != NULL)
                    gN_(i, j, k) = cellVol(i, j + 1, k)/(cellVol(i, j + 1, k) + cellVol(i, j, k));
                else
                    gN_(i, j, k) = 0.;

                if(j > 0 || southBoundaryMeshPtr_ != NULL)
                    gS_(i, j, k) = cellVol(i, j - 1, k)/(cellVol(i, j - 1, k) + cellVol(i, j, k));
                else
                    gS_(i, j, k) = 0.;

                if(k < uCellK() || topBoundaryMeshPtr_ != NULL)
                    gT_(i, j, k) = cellVol(i, j, k + 1)/(cellVol(i, j, k + 1) + cellVol(i, j, k));
                else
                    gT_(i, j, k) = 0.;

                if(k > 0 || bottomBoundaryMeshPtr_ != NULL)
                    gB_(i, j, k) = cellVol(i, j, k - 1)/(cellVol(i, j, k - 1) + cellVol(i, j, k));
                else
                    gB_(i, j, k) = 0.;
            }
        }
    }
}
