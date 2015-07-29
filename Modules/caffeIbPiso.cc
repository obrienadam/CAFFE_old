#include <iostream>

#include "Input.h"
#include "RunControl.h"
#include "IbPiso.h"
#include "InitialConditions.h"
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

        IbPiso piso(input, mesh);

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

    Parallel::finalize();

    return 0;
}
