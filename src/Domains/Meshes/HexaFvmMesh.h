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
    Array3D<Vector3D> faceNormalsI_;
    Array3D<Vector3D> faceNormalsJ_;
    Array3D<Vector3D> faceNormalsK_;
    //Array3D<Vector3D> distanceVectorsI_;
    //Array3D<Vector3D> distanceVectorsJ_;
    //Array3D<Vector3D> distanceVectorsK_;
    Array3D<double> faceAreasI_;
    Array3D<double> faceAreasJ_;
    Array3D<double> faceAreasK_;
    //Array3D<double> cellDistancesI_;
    //Array3D<double> cellDistancesJ_;
    //Array3D<double> cellDistancesK_;

    //- For accessing a field pointer by name

    std::map< std::string, Field<double>* > scalarFieldRegistry_;
    std::map< std::string, Field<Vector3D>* > vectorFieldRegistry_;

public:

    HexaFvmMesh(){}

    //- Fields

    std::vector< Field<double> > scalarFields;
    std::vector< Field<Vector3D> > vectorFields;

    //- Flux fields

    std::vector < Field<double> > scalarFluxFieldsI, scalarFluxFieldsJ, scalarFluxFieldsK;
    std::vector < Field<Vector3D> > vectorFluxFieldsI, vectorFluxFieldsJ, vectorFluxFieldsK;

    //- Initialization

    void initialize(Input &input);

    void addScalarField(std::string scalarFieldName, int type = AUXILLARY);
    void addVectorField(std::string vectorFieldName, int type = AUXILLARY);

    //- Access

    Field<double>* findScalarField(std::string fieldName);
    Field<Vector3D>* findVectorField(std::string fieldName);

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
    /* Vector3D nesE(int i, int j, int k){ return distanceVectorsI_(i + 1, j, k); }
    Vector3D nesW(int i, int j, int k){ return distanceVectorsI_(i, j, k); }
    Vector3D nesN(int i, int j, int k){ return distanceVectorsJ_(i, j + 1, k); }
    Vector3D nesS(int i, int j, int k){ return distanceVectorsJ_(i, j, k); }
    Vector3D nesT(int i, int j, int k){ return distanceVectorsK_(i, j, k + 1); }
    Vector3D nesB(int i, int j, int k){ return distanceVectorsK_(i, j, k); }
    double cellDistanceE(int i, int j, int k){ return cellDistancesI_(i + 1, j, k); }
    double cellDistanceW(int i, int j, int k){ return cellDistancesI_(i, j, k); }
    double cellDistanceN(int i, int j, int k){ return cellDistancesJ_(i, j + 1, k); }
    double cellDistanceS(int i, int j, int k){ return cellDistancesJ_(i, j, k); }
    double cellDistanceT(int i, int j, int k){ return cellDistancesK_(i, j, k + 1); }
    double cellDistanceB(int i, int j, int k){ return cellDistancesK_(i, j, k); } */

    int nCellsI(){ return cellCenters_.sizeI(); }
    int nCellsJ(){ return cellCenters_.sizeJ(); }
    int nCellsK(){ return cellCenters_.sizeK(); }

    //- Dump mesh data to a text file for debugging

    void writeDebug();

};

#endif
