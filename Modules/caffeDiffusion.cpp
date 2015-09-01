#include <iostream>

#include "Parallel.h"
#include "Output.h"
#include "Input.h"
#include "RunControl.h"

#include "ParallelHexaFvmMesh.h"
#include "Diffusion.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Parallel::initialize();

    Input input;
    RunControl runControl;
    ParallelHexaFvmMesh mesh;

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
                mesh.writeTec360(runControl.simTime(), "solution/");
        }
        runControl.displayEndMessage();

        mesh.writeTec360(runControl.simTime(), "solution/");
    }
    catch(const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    Parallel::finalize();

    return 0;
}
