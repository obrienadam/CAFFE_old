#include <iostream>

#include "Output.h"
#include "ArgsList.h"
#include "Input.h"
#include "RunControl.h"

#include "SolverIncludes.h"
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

        ArgsList args(argc, argv);
        Input input(args.inputFilename);
        RunControl runControl;
        Euler solver;
        HexaFvmMesh mesh;
        Diffusion diffusion;
        LinearAdvection linearAdvection;

        mesh.addScalarField("phi", CONSERVED);
        mesh.addScalarField("mu", AUXILLARY);
        mesh.addVectorField("v", AUXILLARY);
        mesh.addScalarField("rho", AUXILLARY);

        // Initialize objects

        runControl.initialize(input);
        solver.initialize(input);
        mesh.initialize(input);
        diffusion.initialize(input, mesh, "phi");
        linearAdvection.initialize(input, mesh, "phi", "v");

        // Set the boundary conditions

        mesh.findScalarField("phi").setAllBoundaries(FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 1.,
                                                     FIXED, 1.,
                                                     FIXED, 1.,
                                                     FIXED, 0.);

        mesh.writeDebug();

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {
            runControl.residualNorm = solver.solve(runControl.timeStep(), diffusion);
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
