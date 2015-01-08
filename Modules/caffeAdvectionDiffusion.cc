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
        LinearAdvection linearAdvection;

        mesh.addScalarField("phi", CONSERVED);
        mesh.addScalarField("mu", AUXILLARY);
        mesh.addVectorField("v", AUXILLARY);
        mesh.addScalarField("rho", AUXILLARY);

        // Initialize objects

        runControl.initializeCase(solver, mesh);
        diffusion.initialize(mesh, "phi");
        linearAdvection.initialize(mesh, "phi", "v");

        // Set boundaries

        mesh.findScalarField("phi").setAllBoundaries(ZERO_GRADIENT, 0.,
                                                     FIXED, 1.,
                                                     FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 0.);

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
