#include <iostream>
#include <string>

#include "HexaMeshGen.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    try
    {

        HexaMeshGen hexaMeshGen(argc, argv);

        hexaMeshGen.readMeshInputFile();

        hexaMeshGen.generateMesh();

        hexaMeshGen.writeMeshFile();

    }
    catch (const char* errorMessage)
    {

        cerr << "Error: " << errorMessage << endl;

    }

    return 0;
    
}
