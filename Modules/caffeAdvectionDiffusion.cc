#include <iostream>

#include "Output.h"
#include "RunControl.h"

#include "Euler.h"
#include "HexaFvmMesh.h"
#include "Diffusion.h"
#include "LinearAdvection.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Output::printCaffeHeader();

    try
    {
        // Declare the basic program objects

        RunControl runControl(argc, argv);
        Euler solver;
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

        // Set the boundary conditions

        mesh.findScalarField("phi").setAllBoundaries(FIXED, 1.,
                                                     FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 0.);

        mesh.writeDebug();

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {
            solver.solve(1e-3, diffusion);
            runControl.displayUpdateMessage();
        }

        // Display end message, the run ended normally

        runControl.displayEndMessage();
        mesh.writeTec360(0.);
    }

    // Catch any exceptions thrown during the run

    catch(const char* errorMessage)
    {
        cerr << "Error: " << errorMessage << endl;
    }

    return 0;
}
