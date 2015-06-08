#include <iostream>

#include "Output.h"
#include "Input.h"
#include "RunControl.h"

#include "SolverIncludes.h"
#include "HexaFvmMesh.h"
#include "Simple.h"
#include "InitialConditions.h"

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
        InitialConditions initialConditions;


        mesh.addVectorField("u", CONSERVED);
        mesh.addScalarField("p", CONSERVED);
        mesh.addScalarField("rho", PRIMITIVE);
        mesh.addScalarField("mu", PRIMITIVE);
        mesh.addScalarField("massFlow", PRIMITIVE);

        // Initialize objects

        runControl.initialize(input);
        solver.initialize(input);
        mesh.initialize(input);
        simple.initialize(input, mesh);
        initialConditions.initialize(mesh);

        initialConditions.readInputFile("case/initialConditions.in");

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {
            runControl.residualNorm = solver.solve(runControl.timeStep(), simple);

            if(runControl.writeToScreen())
            {
                runControl.displayUpdateMessage();
                simple.displayUpdateMessage();
            }

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
