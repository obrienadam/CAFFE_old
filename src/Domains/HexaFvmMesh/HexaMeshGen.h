#ifndef HEXA_MESH_GEN_H
#define HEXA_MESH_GEN_H

#include <vector>
#include <fstream>
#include <string>

#include "Array3D.h"
#include "Point3D.h"

class HexaMeshGen
{

public:

    //- Constructors and destructors
    HexaMeshGen();

    //- HexaMesh file input and output
    void readMeshInputFile(std::string filename);
    void writeMeshFile();

    //- Mesh generation
    void generateMesh();
    void generateBoxMesh(double dx, double dy, double dz);

    //- Mesh check
    void checkMesh();

private:

    //- Helper functions for reading mesh files
    void readVertices(std::ifstream& inputFile);
    void readResolution(std::ifstream& inputFile);

    double metricConversion_;

    std::vector<Point3D> vertices_;
    Array3D<Point3D> nodes_;
    std::string meshName_;
};

#endif
