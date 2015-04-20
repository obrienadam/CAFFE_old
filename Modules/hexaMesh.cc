#include <iostream>
#include <string>

#include "HexaMeshGen.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    try
    {
        HexaMeshGen hexaMeshGen;
        hexaMeshGen.readMeshInputFile("mesh/structuredMesh.in");
        hexaMeshGen.generateMesh();
        hexaMeshGen.writeMeshFile();
    }
    catch (const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    return 0;
}
