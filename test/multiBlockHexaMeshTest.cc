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
        mesh.writeTec360(0., "solution");

        if(Parallel::processNo() == 0)
        {
            cout << mesh().cellXc(11, 4, 4) << endl;
        }
    }
    catch(const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    Parallel::finalize();

    return 0;
}
