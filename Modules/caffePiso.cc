#include <iostream>

#include "Input.h"
#include "RunControl.h"
#include "Piso.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    RunControl runControl;
    HexaFvmMesh mesh;

    try
    {
        input.openInputFile("case/case.in");
        runControl.initialize(input);
        mesh.initialize(input);
        Output::print(mesh.meshStats());

        Piso piso(input, mesh);

        runControl.displayStartMessage();
        while(runControl.continueRun())
        {
            runControl.residualNorm = piso.solve(runControl.timeStep());
            piso.displayUpdateMessage();
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
