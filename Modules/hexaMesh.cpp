#include <iostream>
#include <string>

#include "HexaMeshGen.h"
#include "Parallel.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Parallel::initialize();

    HexaMeshGen hexaMeshGen;

    try
    {
        hexaMeshGen.readFile();
        hexaMeshGen.generateMesh();
        hexaMeshGen.checkMesh();
        hexaMeshGen.writeMeshFile();
    }
    catch (const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    Parallel::finalize();

    return 0;
}
