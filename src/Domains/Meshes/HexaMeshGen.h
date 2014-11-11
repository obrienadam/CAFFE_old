#ifndef HEXA_MESH_GEN_H
#define HEXA_MESH_GEN_H

#include <vector>
#include <fstream>

#include <boost/program_options.hpp>

#include "ArgsList.h"
#include "Array3D.h"
#include "Point3D.h"

class HexaMeshGen
{

private:

    ArgsList argsList_;

    double metricConversion_;

    std::vector<Point3D> vertices_;
    Array3D<Point3D> nodes_;

    void processBuffer(std::string& buffer, bool removeAllWhitespace = true);

    void readVertices(std::ifstream& inputFile);
    void readResolution(std::ifstream& inputFile);
    double getNextElement(std::string& buffer);

public:

    //- Constructors and destructors

    HexaMeshGen();
    HexaMeshGen(int argc, const char* argv[]);

    //- HexaMesh file input and output

    void readMeshInputFile();
    void writeMeshFile();

    //- Mesh generation

    void generateMesh();
    void generateBoxMesh(double dx, double dy, double dz);

};

#endif
