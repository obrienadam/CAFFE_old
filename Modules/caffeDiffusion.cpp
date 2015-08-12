#include <iostream>

#include "Output.h"
#include "Input.h"
#include "RunControl.h"

#include "HexaFvmMesh.h"
#include "Diffusion.h"
#include "Parallel.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    RunControl runControl;
    HexaFvmMesh mesh;

    Parallel::initialize();

    try
    {
        runControl.initialize(input);
        mesh.initialize("mesh/structuredMesh.dat");
        Output::print(mesh.meshStats());

        Diffusion diffusion(input, mesh);

        runControl.displayStartMessage();
        while(runControl.continueRun())
        {
            runControl.residualNorm = diffusion.solve(runControl.timeStep());
            diffusion.displayUpdateMessage();
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

    Parallel::finalize();

    return 0;
}
