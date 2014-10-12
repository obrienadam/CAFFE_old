#include <iostream>

#include "Output.h"
#include "RunControl.h"
#include "SmartPointer.h"

#include "Euler.h"

#include "DomainIncludes.h"
#include "SolverIncludes.h"
#include "SchemeIncludes.h"

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


        //- Check if the command line arguments are valid

        if(!(runControl.validOptionsSelected()))
            return 0;

        SmartPointer<DomainInterface> domain = new HexaFdmMesh;
        SmartPointer<SolverInterface> solver = new Euler;
        SmartPointer<SchemeInterface> scheme = new FiniteDifference;

        //- Initialize the module specific fields

        ScalarField phi("phi");
        VectorField alpha("alpha");
        VectorField a("a");

        domain->addField(phi);;
        domain->addAuxField(alpha);
        domain->addAuxField(a);

        //- Initialize the module objects from the user input

        runControl.initializeObjects(domain,
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
