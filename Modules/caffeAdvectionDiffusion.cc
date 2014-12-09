#include <iostream>

#include "Output.h"
#include "RunControl.h"

#include "Solver.h"
#include "HexaFvmMesh.h"
#include "Diffusion.h"
#include "LinearAdvection.h"

int main(int argc, const char* argv[])
{

    using namespace std;

    Output::displayCaffeHeader();

    try
    {

        // Declare the basic program objects

        RunControl runControl(argc, argv);
        Solver solver;
        HexaFvmMesh mesh;
        Diffusion diffusion;
        LinearAdvection linearadvection;

        mesh.addScalarField("phi", CONSERVED);
        mesh.addVectorField("mu", AUXILLARY);
        mesh.addVectorField("a", AUXILLARY);

        // Initialize objects

        runControl.initializeCase(solver, mesh);
        diffusion.setMeshPointer(&mesh);
        linearadvection.setMeshPointer(&mesh);

        mesh.writeDebug();

        // Display a start message and begin the run

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {

            // runControl.solver.integrateTime();

            // Display a periodic solution update

            runControl.displayUpdateMessage();

        }

        // Display end message, the run ended normally

        runControl.displayEndMessage();

    }

    // Catch any exceptions thrown during the run

    catch(const char* errorMessage)
    {

        cerr << "Error: " << errorMessage << endl;

    }

    return 0;

}
