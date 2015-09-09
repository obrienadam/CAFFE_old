#include <iostream>
#include <memory>

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
    unique_ptr<HexaFvmMesh> meshPtr;

    if(Parallel::nProcesses() == 1)
        meshPtr = unique_ptr<HexaFvmMesh>(new HexaFvmMesh);
    else
        meshPtr = unique_ptr<HexaFvmMesh>(new ParallelHexaFvmMesh);

    try
    {
        runControl.initialize(input);
        meshPtr->initialize("mesh/structuredMesh.dat");
        Output::print(meshPtr->meshStats());

        Diffusion diffusion(input, *meshPtr);

        runControl.displayStartMessage();
        while(runControl.continueRun())
        {
            runControl.residualNorm = diffusion.solve(runControl.timeStep());
            diffusion.displayUpdateMessage();
            runControl.displayUpdateMessage();

            if(runControl.writeToFile())
                meshPtr->writeTec360(runControl.simTime(), "solution/");
        }
        runControl.displayEndMessage();

        meshPtr->writeTec360(runControl.simTime(), "solution/");
    }
    catch(const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    Parallel::finalize();

    return 0;
}
