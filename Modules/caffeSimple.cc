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

        mesh.addVectorField("u", CONSERVED);
        mesh.addScalarField("p", CONSERVED);

        // Initialize objects

        runControl.initialize(input);
        solver.initialize(input);
        mesh.initialize(input);
        simple.initialize(input, mesh);

        // Set the boundary conditions

        mesh.findVectorField("u").setAllBoundaries(ZERO_GRADIENT, Vector3D(0., 0., 0.),
                                                   FIXED, Vector3D(1., 0., 0.),
                                                   FIXED, Vector3D(0., 0., 0.),
                                                   FIXED, Vector3D(0., 0., 0.),
                                                   FIXED, Vector3D(0., 0., 0.),
                                                   FIXED, Vector3D(0., 0., 0.));

        mesh.findScalarField("p").setAllBoundaries(ZERO_GRADIENT, 0.,
                                                   ZERO_GRADIENT, 0.,
                                                   ZERO_GRADIENT, 0.,
                                                   ZERO_GRADIENT, 0.,
                                                   ZERO_GRADIENT, 0.,
                                                   ZERO_GRADIENT, 0.);

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {
            runControl.residualNorm = solver.solve(runControl.timeStep(), simple);

            if(runControl.writeToScreen())
                runControl.displayUpdateMessage();

            if(runControl.writeToFile())
                mesh.writeTec360(runControl.simTime(), "solution");
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
