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
#include <map>

#include "StructuredMesh.h"
#include "Array3D.h"
#include "Vector3D.h"
#include "Point3D.h"
#include "Field.h"

enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};
enum Ordering{ROW, COLUMN, LAYER};

class HexaFvmMesh : public StructuredMesh
{
private:

    //- Geometric data pertaining only to the interior fvm mesh

    Array3D<Point3D> cellCenters_;
    Array3D<double> cellVolumes_;

    Array3D<Point3D> faceCentersI_;
    Array3D<Point3D> faceCentersJ_;
    Array3D<Point3D> faceCentersK_;
    Array3D<double> faceAreasI_;
    Array3D<double> faceAreasJ_;
    Array3D<double> faceAreasK_;
    Array3D<Vector3D> faceNormalsI_;
    Array3D<Vector3D> faceNormalsJ_;
    Array3D<Vector3D> faceNormalsK_;

    Array3D<Vector3D> cellToCellDistanceVectorsI_;
    Array3D<Vector3D> cellToCellDistanceVectorsJ_;
    Array3D<Vector3D> cellToCellDistanceVectorsK_;
    Array3D<double> cellToCellDistancesI_;
    Array3D<double> cellToCellDistancesJ_;
    Array3D<double> cellToCellDistancesK_;

    Array3D<Vector3D> cellToFaceDistanceVectorsE_;
    Array3D<Vector3D> cellToFaceDistanceVectorsW_;
    Array3D<Vector3D> cellToFaceDistanceVectorsN_;
    Array3D<Vector3D> cellToFaceDistanceVectorsS_;
    Array3D<Vector3D> cellToFaceDistanceVectorsT_;
    Array3D<Vector3D> cellToFaceDistanceVectorsB_;
    Array3D<double> cellToFaceDistancesE_;
    Array3D<double> cellToFaceDistancesW_;
    Array3D<double> cellToFaceDistancesN_;
    Array3D<double> cellToFaceDistancesS_;
    Array3D<double> cellToFaceDistancesT_;
    Array3D<double> cellToFaceDistancesB_;

    //- Global cell index map, for implicit methods

    Array3D<int> rowVectorOrdering_;
    Array3D<int> columnVectorOrdering_;
    Array3D<int> layerVectorOrdering_;

    //- Private helper methods

    void initializeCells();
    void initializeCellToCellParameters();
    void initializeFaces();
    void initializeCellToFaceParameters();
    void initializeGlobalIndexMaps();

public:

    HexaFvmMesh(){}

    /**
     * @brief scalarFields A vector containing all scalar fields defined on the domain.
     */
    std::vector< Field<double> > scalarFields;

    /**
     * @brief vectorFields A vector containing all vector fields defined on the domain.
     */
    std::vector< Field<Vector3D> > vectorFields;

    /**
     * @brief initialize Initialize the domain.
     * @param input Input object containg initialization data.
     */
    void initialize(Input &input);

    /**
     * @brief addScalarField Add a new scalar field to the domain.
     * @param scalarFieldName Name of the new scalar field.
     * @param type The type of field, either CONSERVED or AUXILLARY.
     */
    void addScalarField(std::string scalarFieldName, int type = AUXILLARY);

    /**
     * @brief addVectorField Add a new vector field to the domain.
     * @param vectorFieldName Name of the new vector field.
     * @param type The type of field, either CONSERVED or AUXILLARY.
     */
    void addVectorField(std::string vectorFieldName, int type = AUXILLARY);

    /**
     * @brief findScalarField Locate a scalar field within the domain.
     * @param fieldName The name of the field to be searched for.
     * @return A reference to the field if found.
     */
    Field<double>& findScalarField(const std::string &fieldName);

    /**
     * @brief findVectorField Locate a vector field within the domain.
     * @param fieldName The name of the field to be searched for.
     * @return  A reference to  the field if found.
     */
    Field<Vector3D>& findVectorField(const std::string &fieldName);

    Point3D cellXc(int i, int j, int k){ return cellCenters_(i, j, k); }
    double cellVol(int i, int j, int k){ return cellVolumes_(i, j, k); }

    Point3D faceXcE(int i, int j, int k){ return faceCentersI_(i + 1, j, k); }
    Point3D faceXcW(int i, int j, int k){ return faceCentersI_(i, j, k); }
    Point3D faceXcN(int i, int j, int k){ return faceCentersJ_(i, j + 1, k); }
    Point3D faceXcS(int i, int j, int k){ return faceCentersJ_(i, j, k); }
    Point3D faceXcT(int i, int j, int k){ return faceCentersK_(i, j, k + 1); }
    Point3D faceXcB(int i, int j, int k){ return faceCentersK_(i, j, k); }

    double faceAreaE(int i, int j, int k){ return faceAreasI_(i + 1, j, k); }
    double faceAreaW(int i, int j, int k){ return faceAreasI_(i, j, k); }
    double faceAreaN(int i, int j, int k){ return faceAreasJ_(i, j + 1, k); }
    double faceAreaS(int i, int j, int k){ return faceAreasJ_(i, j, k); }
    double faceAreaT(int i, int j, int k){ return faceAreasK_(i, j, k + 1); }
    double faceAreaB(int i, int j, int k){ return faceAreasK_(i, j, k); }

    Vector3D nesE(int i, int j, int k);
    Vector3D nesW(int i, int j, int k);
    Vector3D nesN(int i, int j, int k);
    Vector3D nesS(int i, int j, int k);
    Vector3D nesT(int i, int j, int k);
    Vector3D nesB(int i, int j, int k);

    double cellToCellDistanceE(int i, int j, int k);
    double cellToCellDistanceW(int i, int j, int k);
    double cellToCellDistanceN(int i, int j, int k);
    double cellToCellDistanceS(int i, int j, int k);
    double cellToCellDistanceT(int i, int j, int k);
    double cellToCellDistanceB(int i, int j, int k);

    Vector3D nefE(int i, int j, int k){ return cellToFaceDistanceVectorsE_(i, j, k); }
    Vector3D nefW(int i, int j, int k){ return cellToFaceDistanceVectorsW_(i, j, k); }
    Vector3D nefN(int i, int j, int k){ return cellToFaceDistanceVectorsN_(i, j, k); }
    Vector3D nefS(int i, int j, int k){ return cellToFaceDistanceVectorsS_(i, j, k); }
    Vector3D nefT(int i, int j, int k){ return cellToFaceDistanceVectorsT_(i, j, k); }
    Vector3D nefB(int i, int j, int k){ return cellToFaceDistanceVectorsB_(i, j, k); }

    double cellToFaceDistanceE(int i, int j, int k){ return cellToFaceDistancesE_(i, j, k); }
    double cellToFaceDistanceW(int i, int j, int k){ return cellToFaceDistancesW_(i, j, k); }
    double cellToFaceDistanceN(int i, int j, int k){ return cellToFaceDistancesN_(i, j, k); }
    double cellToFaceDistanceS(int i, int j, int k){ return cellToFaceDistancesS_(i, j, k); }
    double cellToFaceDistanceT(int i, int j, int k){ return cellToFaceDistancesT_(i, j, k); }
    double cellToFaceDistanceB(int i, int j, int k){ return cellToFaceDistancesB_(i, j, k); }

    /**
     * @brief globalIndex Get the vector index for a cell, useful for assembling matrices.
     * @param vectorOrdering The ordering of the vector, either by ROW, COLUMN or LAYER.
     * @return The global index.
     */
    int globalIndex(int i, int j, int k, Ordering vectorOrdering);

    int nCellsI(){ return cellCenters_.sizeI(); }
    int nCellsJ(){ return cellCenters_.sizeJ(); }
    int nCellsK(){ return cellCenters_.sizeK(); }

    /**
     * @brief Output all mesh data to a file for debugging purposes.
     */
    void writeDebug();

    /**
     * @brief Write cell-centered data to the Tecplot360 ASCII format.
     * @param time The solution time.
     */
    void writeTec360(double time = 0);
};

#endif
