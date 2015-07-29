#include <iostream>
#include <string>

#include "MultiBlockHexaMeshGen.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    try
    {
        MultiBlockHexaMeshGen multiBlockHexaMeshGen;

        multiBlockHexaMeshGen.writeMeshFiles();
    }
    catch (const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    return 0;
}
