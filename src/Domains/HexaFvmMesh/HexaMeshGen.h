#ifndef HEXA_MESH_GEN_H
#define HEXA_MESH_GEN_H

#include <vector>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "Array3D.h"
#include "Point3D.h"

class HexaMeshGen
{

public:

    //- Constructors and destructors
    HexaMeshGen();

    //- Mesh generation
    void readFile(const std::string &directory);
    void generateMesh();
    void writeMeshFile(const std::string &directory);

    //- Mesh check
    void checkMesh();

private:

    std::string meshName_;
    double metricConversion_;

    Point3D vertices_[8];
    Vector3D resolution_;
    Array3D<Point3D> nodes_;

    boost::property_tree::ptree meshParameters_;
};

#endif
