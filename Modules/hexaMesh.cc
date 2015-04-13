#include <iostream>
#include <string>

#include "ArgsList.h"
#include "HexaMeshGen.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    try
    {
        ArgsList args(argc, argv);
        HexaMeshGen hexaMeshGen;
        hexaMeshGen.readMeshInputFile(args.inputFilename);
        hexaMeshGen.generateMesh();
        hexaMeshGen.writeMeshFile();
    }
    catch (const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    return 0;
}
