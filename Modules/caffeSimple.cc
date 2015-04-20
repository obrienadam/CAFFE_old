#include <iostream>

#include "Output.h"
#include "Input.h"
#include "RunControl.h"

#include "SolverIncludes.h"
#include "HexaFvmMesh.h"
#include "Simple.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Output::printCaffeHeader();

    try
    {
        // Declare the basic program objects

        Input input("case/case.in");
        RunControl runControl;
        Euler solver;
        HexaFvmMesh mesh;
        Simple simple;

        mesh.addScalarField("phi", CONSERVED);
        mesh.addScalarField("mu", AUXILLARY);
        mesh.addVectorField("v", AUXILLARY);
        mesh.addScalarField("rho", AUXILLARY);

        // Initialize objects

        runControl.initialize(input);
        solver.initialize(input);
        mesh.initialize(input);

        // Set the boundary conditions

        mesh.findScalarField("phi").setAllBoundaries(FIXED, 0.,
                                                     FIXED, 0.,
                                                     FIXED, 1.,
                                                     FIXED, 1.,
                                                     FIXED, 1.,
                                                     FIXED, 0.);

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {
            //runControl.residualNorm = solver.solve(runControl.timeStep(), diffusion);

            if(runControl.writeToScreen())
                runControl.displayUpdateMessage();

            if(runControl.writeToFile())
                mesh.writeTec360(runControl.simTime(), "solution");
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
