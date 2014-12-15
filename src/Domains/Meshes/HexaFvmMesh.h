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

    Array3D<Vector3D> cellToFaceDistanceVectorsI_;
    Array3D<Vector3D> cellToFaceDistanceVectorsJ_;
    Array3D<Vector3D> cellToFaceDistanceVectorsK_;
    Array3D<double> cellToFaceDistancesI_;
    Array3D<double> cellToFaceDistancesJ_;
    Array3D<double> cellToFaceDistancesK_;

    //- For accessing a field pointer by name

    std::map< std::string, Field<double>* > scalarFieldRegistry_;
    std::map< std::string, Field<Vector3D>* > vectorFieldRegistry_;

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

    Field<double>& findScalarField(std::string fieldName);
    Field<Vector3D>& findVectorField(std::string fieldName);

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
    Vector3D nesE(int i, int j, int k){ return cellToCellDistanceVectorsI_(i, j, k); }
    Vector3D nesW(int i, int j, int k){ return cellToCellDistanceVectorsI_(i - 1, j, k); }
    Vector3D nesN(int i, int j, int k){ return cellToCellDistanceVectorsJ_(i, j, k); }
    Vector3D nesS(int i, int j, int k){ return cellToCellDistanceVectorsJ_(i, j - 1, k); }
    Vector3D nesT(int i, int j, int k){ return cellToCellDistanceVectorsK_(i, j, k); }
    Vector3D nesB(int i, int j, int k){ return cellToCellDistanceVectorsK_(i, j, k - 1); }
    double cellDistanceE(int i, int j, int k){ return cellToCellDistancesI_(i, j, k); }
    double cellDistanceW(int i, int j, int k){ return cellToCellDistancesI_(i - 1, j, k); }
    double cellDistanceN(int i, int j, int k){ return cellToCellDistancesJ_(i, j, k); }
    double cellDistanceS(int i, int j, int k){ return cellToCellDistancesJ_(i, j - 1, k); }
    double cellDistanceT(int i, int j, int k){ return cellToCellDistancesK_(i, j, k); }
    double cellDistanceB(int i, int j, int k){ return cellToCellDistancesK_(i, j, k - 1); }

    int nCellsI(){ return cellCenters_.sizeI(); }
    int nCellsJ(){ return cellCenters_.sizeJ(); }
    int nCellsK(){ return cellCenters_.sizeK(); }

    //- Dump mesh data to a text file for debugging

    void writeDebug();

};

#endif
