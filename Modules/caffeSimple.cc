#include <iostream>

#include "Input.h"
#include "RunControl.h"
#include "Simple.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    RunControl runControl;
    HexaFvmMesh mesh;

    try
    {
        runControl.initialize(input);
        mesh.initialize(input);
        Output::print(mesh.meshStats());

        Simple simple(input, mesh);

        runControl.displayStartMessage();
        while(runControl.continueRun())
        {
            runControl.residualNorm = simple.solve(runControl.timeStep());
            simple.displayUpdateMessage();
            runControl.displayUpdateMessage();

            if(runControl.writeToFile())
                mesh.writeTec360(runControl.simTime(), "solution");
        }
        runControl.displayEndMessage();

        mesh.writeTec360(runControl.simTime(), "solution");
    }
    catch(const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    return 0;
}
