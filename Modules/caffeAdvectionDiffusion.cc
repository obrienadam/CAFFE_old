#include <iostream>

#include "Output.h"
#include "RunControl.h"
#include "SmartPointer.h"

#include "AdvectionDiffusion.h"
#include "DomainIncludes.h"
#include "SolverIncludes.h"
#include "SchemeIncludes.h"

typedef HexaFdmMesh<AdvectionDiffusion> Mesh;
typedef SolverInterface<Mesh, AdvectionDiffusion> Solver;

int main(int argc, const char* argv[])
{

    using namespace std;

    Output::displayCaffeHeader();

    try
    {

        //- Declare the basic module objects

        RunControl runControl(argc, argv);
        Solver* solver;

        //- Check if the command line arguments are valid

        if(!(runControl.validOptionsSelected()))
            return 0;

        //- Display a start message and begin the run

        runControl.solverInitialize(solver);

        runControl.displayStartMessage();

        while(runControl.continueRun())
        {

            //- Display a periodic solution update

            runControl.displayUpdateMessage();

        }

        //- Display end message, the run ended normally

        runControl.displayEndMessage();

    }

    //- Catch any exceptions thrown during the run

    catch(const char* errorMessage)
    {

        cerr << "Error: " << errorMessage << endl;

    }

    return 0;

}
