#ifndef PARALLEL_H
#define PARALLEL_H

#include <vector>

#include "Vector3D.h"
#include "Field.h"
#include "Input.h"
#include "Matrix.h"

class Parallel
{

private:

    static int commBufferSize_;
    static std::vector<double> commBuffer_;

public:

    static void initialize();
    static void finalize();

    static int nProcesses();
    static int processNo();
    static bool isThisProcessor(int source);
    static bool isMainProcessor();
    static void ownerShipRange(int nEntities, int &iLower, int &iUpper);

    static int sum(int number);
    static int min(int number);
    static int max(int number);

    static double sum(double number);
    static double min(double number);
    static double max(double number);

    static void broadcastInput(int source, Input& input);

    //- Note: All of these methods assume that the data structures have been properly sized!
    static void send(int source, int dest, std::vector<double>& vector);
    static void send(int source, int dest, Array3D<double>& doubleArray3D);
    static void send(int source, int dest, Array3D<Vector3D>& vector3DArray3D);
    static void send(int source, int dest, Matrix& matrix);

    static void barrier();
};

#endif
