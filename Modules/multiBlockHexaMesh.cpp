#include <iostream>
#include <string>

#include "MultiBlockHexaMeshGen.h"
#include "Parallel.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Parallel::initialize();

    MultiBlockHexaMeshGen multiBlockHexaMeshGen;

    try
    {
        multiBlockHexaMeshGen.readFile();
        multiBlockHexaMeshGen.writeMeshFiles();
    }
    catch (const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    Parallel::finalize();

    return 0;
}
