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
    Array3D<double> faceAreasI_;
    Array3D<double> faceAreasJ_;
    Array3D<double> faceAreasK_;

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

    Point3D cellXc(int i, int j, int k);
    double cellVol(int i, int j, int k);

    Point3D faceXcE(int i, int j, int k);
    Point3D faceXcW(int i, int j, int k);
    Point3D faceXcN(int i, int j, int k);
    Point3D faceXcS(int i, int j, int k);
    Point3D faceXcT(int i, int j, int k);
    Point3D faceXcB(int i, int j, int k);
    double faceAreaE(int i, int j, int k);
    double faceAreaW(int i, int j, int k);
    double faceAreaN(int i, int j, int k);
    double faceAreaS(int i, int j, int k);
    double faceAreaT(int i, int j, int k);
    double faceAreaB(int i, int j, int k);

    //- Dump mesh data to a text file for debugging

    void writeDebug();

};

#endif
