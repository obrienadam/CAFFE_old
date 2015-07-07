#include <iostream>
#include <cstdlib>

#include "Parallel.h"
#include "Array3D.h"
#include "Vector3D.h"

int main()
{
    using namespace std;

    int i, j, k;
    Array3D<double> doubleArray3D;
    Array3D<Vector3D> vector3DArray3D;

    Parallel::initialize();

    try
    {

        if(Parallel::isMainProcessor())
            cout << "Running test with " << Parallel::nProcesses() << " processes." << endl;

        Parallel::barrier();

        cout << "Hello world from process number " << Parallel::processNo() << "!" << endl;

        Parallel::barrier();

        if(Parallel::isMainProcessor())
            cout << "Initializing a random double 3D array on process " << Parallel::processNo() << "." << endl;

        doubleArray3D.allocate(100, 100, 100);

        if(Parallel::isMainProcessor())
        {
            for(i = 0; i < doubleArray3D.size(); ++i)
                doubleArray3D(i) = double(rand())/double(RAND_MAX);
        }

        if(Parallel::isMainProcessor())
            cout << "Sending random double 3D array to process 1." << endl;

        Parallel::send(0, 1, doubleArray3D);

        i = 40;
        j = 40;
        k = 40;

        if(Parallel::processNo() == 0)
        {
            cout << "On process " << Parallel::processNo() << ", element " << i << ", " << j << ", " << k << " contains " << doubleArray3D(i, j, k) << endl;
        }
        if(Parallel::processNo() == 1)
        {
            cout << "On process " << Parallel::processNo() << ", element " << i << ", " << j << ", " << k << " contains " << doubleArray3D(i, j, k) << endl;
        }

        Parallel::barrier();

        if(Parallel::isMainProcessor())
            cout << "Initializing a random vector3D 3D array on process " << Parallel::processNo() << "." << endl;

        vector3DArray3D.allocate(100, 100, 100);
        cout << vector3DArray3D(i, j, k) << endl;

        if(Parallel::isMainProcessor())
        {
            for(i = 0; i < vector3DArray3D.size(); ++i)
                vector3DArray3D(i) = Vector3D(double(rand())/double(RAND_MAX), double(rand())/double(RAND_MAX), double(rand())/double(RAND_MAX));
        }

        i = 29;
        j = 39;
        k = 11;

        Parallel::send(0, 1, vector3DArray3D);

        if(Parallel::processNo() == 0)
        {
            cout << "On process " << Parallel::processNo() << ", element " << i << ", " << j << ", " << k << " contains " << vector3DArray3D(i, j, k) << endl;
        }
        if(Parallel::processNo() == 1)
        {
            cout << "On process " << Parallel::processNo() << ", element " << i << ", " << j << ", " << k << " contains " << vector3DArray3D(i, j, k) << endl;
        }

        Parallel::send(1, 2, vector3DArray3D);

        if(Parallel::processNo() == 1)
        {
            cout << "On process " << Parallel::processNo() << ", element " << i << ", " << j << ", " << k << " contains " << vector3DArray3D(i, j, k) << endl;
        }
        if(Parallel::processNo() == 2)
        {
            cout << "On process " << Parallel::processNo() << ", element " << i << ", " << j << ", " << k << " contains " << vector3DArray3D(i, j, k) << endl;
        }

        int nEntities = 9312, iLower, iUpper;

        Parallel::barrier();

        if(Parallel::isMainProcessor())
            cout << endl << "Testing the ownership range calculation with " << nEntities << " entities." << endl;

        Parallel::ownerShipRange(nEntities, iLower, iUpper);

        for(i = 0; i < Parallel::nProcesses(); ++i)
        {
            Parallel::barrier();

            if(Parallel::isThisProcessor(i))
            {
                cout << "Ownership range for processor " << Parallel::processNo() << ": " << iLower << " -- " << iUpper << endl;
            }
        }
    }
    catch(const char* errorMessage)
    {
        cerr << errorMessage << endl;
    }

    Parallel::finalize();
}
