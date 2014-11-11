#ifndef HEXA_MESH_H
#define HEXA_MESH_H

#include <boost/program_options.hpp>

#include "ArgsList.h"
#include "Array3D.h"
#include "Point3D.h"

class HexaMeshGen
{

private:

    ArgsList argsList_;

    Array3D<Point3D> vertices_;
    Array3D<Point3D> nodes_;

public:

    HexaMeshGen();
    HexaMeshGen(int argc, const char* argv[]);

    void readMeshInputFile();
    void writeMeshFile();

    void generateMesh();
    void generateBoxMesh(double dx, double dy, double dz);

};

#endif
