#include <iostream>

#include "Output.h"
#include "Input.h"
#include "RunControl.h"

#include "SolverIncludes.h"
#include "HexaFvmMesh.h"
#include "IbSimple.h"
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
        IbSimple ibSimple;
        InitialConditions initialConditions;

        mesh.addVectorField("u", CONSERVED);
        mesh.addScalarField("p", CONSERVED);
        mesh.addScalarField("rho", PRIMITIVE);
        mesh.addScalarField("mu", PRIMITIVE);
        mesh.addScalarField("massFlow", PRIMITIVE);
        mesh.addScalarField("ibField", AUXILLARY);

        // Initialize objects

        runControl.initialize(input);
        solver.initialize(input);
        mesh.initialize(input);
        ibSimple.initialize(input, mesh);
        initialConditions.initialize(mesh);

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {
            runControl.residualNorm = solver.solve(runControl.timeStep(), ibSimple);

            if(runControl.writeToScreen())
            {
                runControl.displayUpdateMessage();
                ibSimple.displayUpdateMessage();
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
