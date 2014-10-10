#include <iostream>

#include "Output.h"
#include "RunControl.h"
#include "HexaFdmMesh.h"

#include "Euler.h"

#include "FiniteDifference.h"

#include "ScalarField.h"
#include "VectorField.h"

int main(int argc, const char* argv[])
{

    using namespace std;

    Output::displayCaffeHeader();

    try
    {

        //- Declare the basic module objects

        RunControl runControl(argc, argv);
        HexaFdmMesh mesh;
        SolverInterface* solver = new Euler;
        SchemeInterface* scheme = new FiniteDifference;

        //- Initialize the module specific fields

        ScalarField phi("phi");
        ScalarField alpha("alpha");
        VectorField a("a");

        mesh.addField(phi);
        mesh.addField(alpha);
        mesh.addField(a);

        runControl.initializeObjects(&mesh,
                                     solver,
                                     scheme);

        //- Display a start message and begin the run

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
