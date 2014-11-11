#include <fstream>

#include "HexaMeshGen.h"

HexaMeshGen::HexaMeshGen()
    :
      vertices_(2, 2, 2)
{

}

HexaMeshGen::HexaMeshGen(int argc, const char *argv[])
    :
      HexaMeshGen()
{

    argsList_.readArgs(argc, argv);

}

void HexaMeshGen::readMeshInputFile()
{
    using namespace std;

    ifstream inputFile(argsList_.varsMap_["file"].as<string>().c_str());

    if(!inputFile.is_open())
    {

        throw "Mesh input file not found.";

    }

}

void HexaMeshGen::writeMeshFile()
{

}

void HexaMeshGen::generateMesh()
{

}

void HexaMeshGen::generateBoxMesh(double dx, double dy, double dz)
{



}
