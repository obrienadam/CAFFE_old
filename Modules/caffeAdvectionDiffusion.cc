#include <iostream>

#include "Output.h"
#include "RunControl.h"

#include "HexaFvmMesh.h"
#include "Solver.h"
#include "SchemeIncludes.h"

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

        mesh.addScalarField("phi");
        mesh.addVectorField("mu");
        mesh.addVectorField("a");

        // Initialize objects

        runControl.initializeCase(solver, mesh);

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
