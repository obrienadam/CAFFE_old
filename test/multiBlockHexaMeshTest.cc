#include <iostream>

#include "MultiBlockHexaFvmMesh.h"
#include "Parallel.h"

int main()
{
    using namespace std;

    MultiBlockHexaFvmMesh mesh;

    Parallel::initialize();

    try
    {
        mesh.initialize();
    }
    catch(const char* errorMessage)
    {
        cerr << errorMessage << endl;
    }

    Parallel::finalize();

    return 0;
}
