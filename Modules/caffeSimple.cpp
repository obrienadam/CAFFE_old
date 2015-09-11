#include <iostream>

#include "Parallel.h"
#include "Input.h"
#include "RunControl.h"
#include "Simple.h"
#include "ParallelHexaFvmMesh.h"

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

        Simple simple(input, *meshPtr);

        runControl.displayStartMessage();
        while(runControl.continueRun())
        {
            runControl.residualNorm = simple.solve(runControl.timeStep());
            simple.displayUpdateMessage();
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
