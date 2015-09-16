#ifndef PARALLEL_H
#define PARALLEL_H

#include <vector>
#include <mpi.h>

#include "Vector3D.h"
#include "Array3D.h"
#include "ParallelBuffers.h"

class Parallel
{
public:

    static const int PROC_NULL;

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
    static std::string broadcast(const std::string &str, int source);

    static void broadcast(int source, std::vector<double> &doubles);
    static void broadcast(int source, std::vector<Vector3D> &vecs);

    static int sum(int number);
    static double sum(double number);

    static int min(int number);
    static double min(double number);

    static int max(int number);
    static double max(double number);

    //- These methods assume the data structures are already apropriately sized
    static void send(int source, int dest, std::vector<double> &doubles);
    static void send(int source, int dest, std::vector<Vector3D> &vecs);

    //- These methods will also resize the receiving data structure if necessary. Involves an extra communication
    static void send(int source, int dest, const Array3D<int> &sourceArray3D, Array3D<int> &destArray3D);

    //- Non-blocking comms. Must remember to call waitAll to ensure comms are in sync
    static void iSend(int source, int dest, int tag, const std::vector<int> &ints);
    static void iRecv(int source, int dest, int tag, std::vector<int> &ints);

    static void iSend(int source, int dest, int tag, const std::vector<double> &doubles);
    static void iRecv(int source, int dest, int tag, std::vector<double> &doubles);

    static void iSend(int source, int dest, int tag, const std::vector<Vector3D> &vecs);
    static void iRecv(int source, int dest, int tag, std::vector<Vector3D> &vecs);

    static void allGather(int number, std::vector<int> &vector);
    static void allGather(double number, std::vector<double> &vector);

    static void barrier();

    //- Important method that must be called at the end of a cycle of non-blocking communication
    static void waitAll();

private:

    static ParallelBuffers buffers_;

    static std::vector< MPI::Request > requests_;
    static int nProcesses_, procNo_;
    static bool isInitialized_;
};

#endif
