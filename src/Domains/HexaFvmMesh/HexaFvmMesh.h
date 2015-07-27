/**
 * @file    HexaFvmMesh.h
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
 * This file contains the interface for class HexaFvmMesh, which
 * is a domain for solving finite-volume problems on arbitrarily
 * shaped hexahedral cells with planar faces.
 */

#ifndef HEXA_FVM_MESH_H
#define HEXA_FVM_MESH_H

#include <string>
#include <vector>

#include "StructuredMesh.h"
#include "Array3D.h"
#include "Vector3D.h"
#include "Point3D.h"

class HexaFvmMesh : public StructuredMesh
{

public:

    enum Direction{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};
    enum {BSW = 0, BSE = 1, BNE = 2, BNW = 3, TSW = 4, TSE = 5, TNE = 6, TNW = 7};

    HexaFvmMesh();
    HexaFvmMesh(const HexaFvmMesh& other);

    void initialize(Input &input);
    void initialize(Array3D<Point3D>& nodes);

    void addBoundaryMesh(const HexaFvmMesh &mesh, Direction relativeLocation);
    bool eastMeshExists() const;
    bool westMeshExists() const;
    bool northMeshExists() const;
    bool southMeshExists() const;
    bool topMeshExists() const;
    bool bottomMeshExists() const;

    virtual std::string meshStats();

    // Retrieve cell parameters
    Point3D cellXc(int i, int j, int k) const;
    double cellVol(int i, int j, int k) const;

    Point3D node(int i, int j, int k, int nodeNo) const;

    // Retrive face parameters
    Point3D faceXcE(int i, int j, int k) const { return faceCentersI_(i + 1, j, k); }
    Point3D faceXcW(int i, int j, int k) const { return faceCentersI_(i, j, k); }
    Point3D faceXcN(int i, int j, int k) const { return faceCentersJ_(i, j + 1, k); }
    Point3D faceXcS(int i, int j, int k) const { return faceCentersJ_(i, j, k); }
    Point3D faceXcT(int i, int j, int k) const { return faceCentersK_(i, j, k + 1); }
    Point3D faceXcB(int i, int j, int k) const { return faceCentersK_(i, j, k); }

    Vector3D fAreaNormE(int i, int j, int k) const { return faceNormalsI_(i + 1, j, k); }
    Vector3D fAreaNormW(int i, int j, int k) const { return -faceNormalsI_(i, j, k); }
    Vector3D fAreaNormN(int i, int j, int k) const { return faceNormalsJ_(i, j + 1, k); }
    Vector3D fAreaNormS(int i, int j, int k) const { return -faceNormalsJ_(i, j, k); }
    Vector3D fAreaNormT(int i, int j, int k) const { return faceNormalsK_(i, j, k + 1); }
    Vector3D fAreaNormB(int i, int j, int k) const { return -faceNormalsK_(i, j, k); }

    double faceAreaE(int i, int j, int k) const { return faceAreasI_(i + 1, j, k); }
    double faceAreaW(int i, int j, int k) const { return faceAreasI_(i, j, k); }
    double faceAreaN(int i, int j, int k) const { return faceAreasJ_(i, j + 1, k); }
    double faceAreaS(int i, int j, int k) const { return faceAreasJ_(i, j, k); }
    double faceAreaT(int i, int j, int k) const { return faceAreasK_(i, j, k + 1); }
    double faceAreaB(int i, int j, int k) const { return faceAreasK_(i, j, k); }

    Point3D faceXcI(int i, int j, int k) const { return faceCentersI_(i, j, k); }
    Point3D faceXcJ(int i, int j, int k) const { return faceCentersJ_(i, j, k); }
    Point3D faceXcK(int i, int j, int k) const { return faceCentersK_(i, j, k); }

    Vector3D fAreaNormI(int i, int j, int k) const { return faceNormalsI_(i, j, k); }
    Vector3D fAreaNormJ(int i, int j, int k) const { return faceNormalsJ_(i, j, k); }
    Vector3D fAreaNormK(int i, int j, int k) const { return faceNormalsK_(i, j, k); }

    double fAreaI(int i, int j, int k) const { return faceAreasI_(i, j, k); }
    double fAreaJ(int i, int j, int k) const { return faceAreasJ_(i, j, k); }
    double fAreaK(int i, int j, int k) const { return faceAreasK_(i, j, k); }

    // Retrieve cell to cell parameters
    Vector3D rCellE(int i, int j, int k) const;
    Vector3D rCellW(int i, int j, int k) const;
    Vector3D rCellN(int i, int j, int k) const;
    Vector3D rCellS(int i, int j, int k) const;
    Vector3D rCellT(int i, int j, int k) const;
    Vector3D rCellB(int i, int j, int k) const;

    Vector3D rnCellE(int i, int j, int k) const;
    Vector3D rnCellW(int i, int j, int k) const;
    Vector3D rnCellN(int i, int j, int k) const;
    Vector3D rnCellS(int i, int j, int k) const;
    Vector3D rnCellT(int i, int j, int k) const;
    Vector3D rnCellB(int i, int j, int k) const;

    double rCellMagE(int i, int j, int k) const;
    double rCellMagW(int i, int j, int k) const;
    double rCellMagN(int i, int j, int k) const;
    double rCellMagS(int i, int j, int k) const;
    double rCellMagT(int i, int j, int k) const;
    double rCellMagB(int i, int j, int k) const;

    // Retrieve cell to face parameters
    Vector3D rFaceE(int i, int j, int k) const { return cellToFaceRelativeVectorsE_(i, j, k); }
    Vector3D rFaceW(int i, int j, int k) const { return cellToFaceRelativeVectorsW_(i, j, k); }
    Vector3D rFaceN(int i, int j, int k) const { return cellToFaceRelativeVectorsN_(i, j, k); }
    Vector3D rFaceS(int i, int j, int k) const { return cellToFaceRelativeVectorsS_(i, j, k); }
    Vector3D rFaceT(int i, int j, int k) const { return cellToFaceRelativeVectorsT_(i, j, k); }
    Vector3D rFaceB(int i, int j, int k) const { return cellToFaceRelativeVectorsB_(i, j, k); }

    Vector3D rnFaceE(int i, int j, int k) const { return cellToFaceUnitVectorsE_(i, j, k); }
    Vector3D rnFaceW(int i, int j, int k) const { return cellToFaceUnitVectorsW_(i, j, k); }
    Vector3D rnFaceN(int i, int j, int k) const { return cellToFaceUnitVectorsN_(i, j, k); }
    Vector3D rnFaceS(int i, int j, int k) const { return cellToFaceUnitVectorsS_(i, j, k); }
    Vector3D rnFaceT(int i, int j, int k) const { return cellToFaceUnitVectorsT_(i, j, k); }
    Vector3D rnFaceB(int i, int j, int k) const { return cellToFaceUnitVectorsB_(i, j, k); }

    double rFaceMagE(int i, int j, int k) const { return cellToFaceDistancesE_(i, j, k); }
    double rFaceMagW(int i, int j, int k) const { return cellToFaceDistancesW_(i, j, k); }
    double rFaceMagN(int i, int j, int k) const { return cellToFaceDistancesN_(i, j, k); }
    double rFaceMagS(int i, int j, int k) const { return cellToFaceDistancesS_(i, j, k); }
    double rFaceMagT(int i, int j, int k) const { return cellToFaceDistancesT_(i, j, k); }
    double rFaceMagB(int i, int j, int k) const { return cellToFaceDistancesB_(i, j, k); }

    int nCellsI() const { return cellCenters_.sizeI(); }
    int nCellsJ() const { return cellCenters_.sizeJ(); }
    int nCellsK() const { return cellCenters_.sizeK(); }
    int nCells() const { return cellCenters_.size(); }
    int nFacesI() const { return faceCentersI_.sizeI(); }
    int nFacesJ() const { return faceCentersJ_.sizeJ(); }
    int nFacesK() const { return faceCentersK_.sizeK(); }

    int uCellI() const { return cellCenters_.sizeI() - 1; }
    int uCellJ() const { return cellCenters_.sizeJ() - 1; }
    int uCellK() const { return cellCenters_.sizeK() - 1; }

    double dE(int i, int j, int k) const { return dE_(i, j, k); }
    double dW(int i, int j, int k) const { return dW_(i, j, k); }
    double dN(int i, int j, int k) const { return dN_(i, j, k); }
    double dS(int i, int j, int k) const { return dS_(i, j, k); }
    double dT(int i, int j, int k) const { return dT_(i, j, k); }
    double dB(int i, int j, int k) const { return dB_(i, j, k); }

    Vector3D cE(int i, int j, int k) const { return cE_(i, j, k); }
    Vector3D cW(int i, int j, int k) const { return cW_(i, j, k); }
    Vector3D cN(int i, int j, int k) const { return cN_(i, j, k); }
    Vector3D cS(int i, int j, int k) const { return cS_(i, j, k); }
    Vector3D cT(int i, int j, int k) const { return cT_(i, j, k); }
    Vector3D cB(int i, int j, int k) const { return cB_(i, j, k); }

    double gE(int i, int j, int k) const { return gE_(i, j, k); }
    double gW(int i, int j, int k) const { return gW_(i, j, k); }
    double gN(int i, int j, int k) const { return gN_(i, j, k); }
    double gS(int i, int j, int k) const { return gS_(i, j, k); }
    double gT(int i, int j, int k) const { return gT_(i, j, k); }
    double gB(int i, int j, int k) const { return gB_(i, j, k); }

    /**
     * @brief Locate the cell whose center is closest to a point, and store the indices.
     * @param point The point to be considered.
     * @param ii The i index of the cell closest to point.
     * @param jj The j index of the cell closest to point.
     * @param kk The k index of the cell closest to point.
     */
    void locateCell(const Point3D& point, int& ii, int& jj, int& kk) const;

    /**
     * @brief Locate the 8 cells that form a hexahedron with their centers and enclose the point.
     * @param point The point of interest.
     * @param ii Array of the i indices enclosing the point.
     * @param jj Array of the j indices enclosing the point.
     * @param kk Array of the k indices enclosing the point.
     */
    void locateEnclosingCells(const Point3D& point, int ii[], int jj[], int kk[]) const;

    /**
     * @brief Output all mesh data to a file for debugging purposes.
     */
    void writeDebug() const;

    void addArray3DToTecplotOutput(std::string name, const Array3D<double> *array3DPtr) const;
    void addArray3DToTecplotOutput(std::string name, const Array3D<Vector3D> *array3DPtr) const;

    /**
     * @brief Write cell-centered data to the Tecplot360 ASCII format.
     * @param time The solution time.
     */
    void writeTec360(double time = 0, std::string directoryName = "");

private:

    //- Private helper methods
    void initializeCellsAndFaces();
    void initializeCells();
    void initializeCellToCellParameters();
    void initializeFaces();
    void initializeCellToFaceParameters();
    void computeMeshMetrics();

    //- Geometric data pertaining only to the interior fvm mesh
    Array3D<Point3D> cellCenters_;
    Array3D<double> cellVolumes_;

    Array3D<Point3D> faceCentersI_;
    Array3D<Point3D> faceCentersJ_;
    Array3D<Point3D> faceCentersK_;
    Array3D<Vector3D> faceNormalsI_;
    Array3D<Vector3D> faceNormalsJ_;
    Array3D<Vector3D> faceNormalsK_;
    Array3D<double> faceAreasI_;
    Array3D<double> faceAreasJ_;
    Array3D<double> faceAreasK_;
    Array3D<Vector3D> faceUnitNormalsI_;
    Array3D<Vector3D> faceUnitNormalsJ_;
    Array3D<Vector3D> faceUnitNormalsK_;

    Array3D<Vector3D> cellToCellRelativeVectorsI_;
    Array3D<Vector3D> cellToCellRelativeVectorsJ_;
    Array3D<Vector3D> cellToCellRelativeVectorsK_;
    Array3D<Vector3D> cellToCellUnitVectorsI_;
    Array3D<Vector3D> cellToCellUnitVectorsJ_;
    Array3D<Vector3D> cellToCellUnitVectorsK_;
    Array3D<double> cellToCellDistancesI_;
    Array3D<double> cellToCellDistancesJ_;
    Array3D<double> cellToCellDistancesK_;

    Array3D<Vector3D> cellToFaceRelativeVectorsE_;
    Array3D<Vector3D> cellToFaceRelativeVectorsW_;
    Array3D<Vector3D> cellToFaceRelativeVectorsN_;
    Array3D<Vector3D> cellToFaceRelativeVectorsS_;
    Array3D<Vector3D> cellToFaceRelativeVectorsT_;
    Array3D<Vector3D> cellToFaceRelativeVectorsB_;
    Array3D<Vector3D> cellToFaceUnitVectorsE_;
    Array3D<Vector3D> cellToFaceUnitVectorsW_;
    Array3D<Vector3D> cellToFaceUnitVectorsN_;
    Array3D<Vector3D> cellToFaceUnitVectorsS_;
    Array3D<Vector3D> cellToFaceUnitVectorsT_;
    Array3D<Vector3D> cellToFaceUnitVectorsB_;
    Array3D<double> cellToFaceDistancesE_;
    Array3D<double> cellToFaceDistancesW_;
    Array3D<double> cellToFaceDistancesN_;
    Array3D<double> cellToFaceDistancesS_;
    Array3D<double> cellToFaceDistancesT_;
    Array3D<double> cellToFaceDistancesB_;

    //- Misc mesh metrics commonly used in fvm
    Array3D<double> dE_, dW_, dN_, dS_, dT_, dB_;
    Array3D<Vector3D> cE_, cW_, cN_, cS_, cT_, cB_;

    //- Volume based interpolation factors
    Array3D<double> gE_, gW_, gN_, gS_, gT_, gB_;

    const HexaFvmMesh *eastBoundaryMeshPtr_;
    const HexaFvmMesh *westBoundaryMeshPtr_;
    const HexaFvmMesh *northBoundaryMeshPtr_;
    const HexaFvmMesh *southBoundaryMeshPtr_;
    const HexaFvmMesh *topBoundaryMeshPtr_;
    const HexaFvmMesh *bottomBoundaryMeshPtr_;

    mutable std::vector< std::pair< std::string, const Array3D<double> *> > scalarVariablePtrs_;
    mutable std::vector< std::pair< std::string, const Array3D<Vector3D> *> > vectorVariablePtrs_;
};

#endif
