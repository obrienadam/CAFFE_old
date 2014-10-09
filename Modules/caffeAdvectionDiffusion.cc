#include <iostream>

#include "RunControl.h"
#include "HexaFdmMesh.h"

#include "Euler.h"

#include "FiniteDifference.h"

#include "ScalarField.h"
#include "VectorField.h"

int main(int argc, const char* argv[])
{

    using namespace std;

    try
    {

        //- Declare the basic module objects

        RunControl runControl(argc, argv);
        HexaFdmMesh mesh(30, 30, 30);
        SolverInterface* solver = new Euler;
        SchemeInterface* scheme = new FiniteDifference;

        //- Initialize the module specific fields

        ScalarField phi("phi", 30, 30, 30);
        ScalarField alpha("alpha", 30, 30, 30);
        VectorField a("a", 30, 30, 30);

        mesh.addField(phi);
        mesh.addField(alpha);
        mesh.addField(a);

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
