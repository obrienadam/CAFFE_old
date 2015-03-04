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

    //- Fields

    std::vector< Field<double> > scalarFields;
    std::vector< Field<Vector3D> > vectorFields;

    //- Initialization

    void initialize(Input &input);

    void addScalarField(std::string scalarFieldName, int type = AUXILLARY);
    void addVectorField(std::string vectorFieldName, int type = AUXILLARY);

    //- Access

    Field<double>& findScalarField(const std::string &fieldName);
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

    int globalIndex(int i, int j, int k, Ordering vectorOrdering);

    int nCellsI(){ return cellCenters_.sizeI(); }
    int nCellsJ(){ return cellCenters_.sizeJ(); }
    int nCellsK(){ return cellCenters_.sizeK(); }

    //- Dump mesh data to a text file for debugging

    void writeDebug();
};

#endif
