#ifndef HEXA_FVM_MESH_H
#define HEXA_FVM_MESH_H

#include <string>
#include <vector>

#include "StructuredMesh.h"
#include "Array3D.h"
#include "Vector3D.h"
#include "Point3D.h"
#include "Field.h"

class HexaFvmMesh : public StructuredMesh
{

private:

    //- Geometric data pertaining only to the fvm mesh

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

public:

    //- Fields

    std::vector< Field<double> > scalarFields_;
    std::vector< Field<Vector3D> > vectorFields_;

public:

    HexaFvmMesh();

    //- Initialization

    void initialize(Input &input);

    void addScalarField(std::string scalarFieldName);
    void addVectorField(std::string vectorFieldName);

    //- Dump mesh data to a text file for debugging

    void writeDebug();

};

#endif
