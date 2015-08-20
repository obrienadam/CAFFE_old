#ifndef PARALLEL_H
#define PARALLEL_H

#include <vector>

#include "Vector3D.h"
#include "Field.h"
#include "Input.h"
#include "Matrix.h"

class Parallel
{
public:

    static void initialize();
    static void finalize();

    static int nProcesses();
    static int processNo();
    static int mainProcNo(){ return 0; }
    static bool isMainProcessor();
    static std::pair<int, int> ownershipRange(int nEntities);

    static int broadcast(int number, int source);
    static double broadcast(double number, int source);
    static Vector3D broadcast(const Vector3D &vec, int source);

    static void send(int source, int dest, std::vector<double> &doubles);
    static void send(int source, int dest, std::vector<Vector3D> &vecs);

    static void allGather(int number, std::vector<int> &vector);
    static void allGather(double number, std::vector<double> &vector);

    static void barrier();

private:

    static void loadBuffer(const std::vector<Vector3D> &vecs);
    static void unloadBuffer(std::vector<Vector3D> &vecs);

    static std::vector<double> commBuffer_;
    static int nProcesses_, procNo_;
    static bool isInitialized_;
};

#endif
