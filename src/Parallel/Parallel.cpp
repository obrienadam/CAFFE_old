#include <mpi.h>

#include "Parallel.h"

/********************** Public methods ***********************/

void Parallel::initialize()
{
    MPI::Init();
    procNo_ = MPI::COMM_WORLD.Get_rank();
    nProcesses_ = MPI::COMM_WORLD.Get_size();
    commBuffer_.resize(100000);
    isInitialized_ = true;
}

void Parallel::finalize()
{
    MPI::Finalize();
    procNo_ = 0;
    nProcesses_ = 1;
    commBuffer_.erase(commBuffer_.begin(), commBuffer_.end());
    isInitialized_ = false;
}

int Parallel::nProcesses()
{
    return nProcesses_;
}

int Parallel::processNo()
{
    return procNo_;
}

bool Parallel::isMainProcessor()
{
    return procNo_ == 0;
}

std::pair<int, int> Parallel::ownershipRange(int nEntities)
{
    int nRemainingEntities, nEntitiesThisProc;
    std::pair<int, int> range;

    nEntitiesThisProc = nEntities/nProcesses();
    nRemainingEntities = nEntities%nProcesses();

    if(processNo() < nRemainingEntities)
    {
        ++nEntitiesThisProc;
        range.first = processNo()*nEntitiesThisProc;
    }
    else
    {
        range.first = (nEntitiesThisProc + 1)*nRemainingEntities + (processNo() - nRemainingEntities)*nEntitiesThisProc;
    }

    range.second = range.first + nEntitiesThisProc - 1;

    return range;
}

int Parallel::broadcast(int number, int source)
{
    MPI::COMM_WORLD.Bcast(&number, 1, MPI::INT, source);
    return number;
}

double Parallel::broadcast(double number, int source)
{
    MPI::COMM_WORLD.Bcast(&number, 1, MPI::DOUBLE, source);
    return number;
}

Vector3D Parallel::broadcast(const Vector3D &vec, int source)
{
    double buff[3] = {vec.x, vec.y, vec.z};
    MPI::COMM_WORLD.Bcast(buff, 3, MPI::DOUBLE, source);
    return Vector3D(buff[0], buff[1], buff[2]);
}

void Parallel::send(int source, int dest, std::vector<double> &doubles)
{
    if(processNo() == source)
    {
        MPI::COMM_WORLD.Send(doubles.data(), doubles.size(), MPI::DOUBLE, dest, 0);
    }
    else if(processNo() == dest)
    {
        MPI::COMM_WORLD.Recv(doubles.data(), doubles.size(), MPI::DOUBLE, source, 0);
    }
}

void Parallel::send(int source, int dest, std::vector<Vector3D> &vecs)
{
    if(processNo() == source)
    {
        loadBuffer(vecs);
        MPI::COMM_WORLD.Send(commBuffer_.data(), 3*vecs.size(), MPI::DOUBLE, dest, 0);
    }
    else if(processNo() == dest)
    {
        MPI::COMM_WORLD.Recv(commBuffer_.data(), 3*vecs.size(), MPI::DOUBLE, source, 0);
        unloadBuffer(vecs);
    }
}

void Parallel::allGather(int number, std::vector<int> &vector)
{
    if(isInitialized_)
        MPI::COMM_WORLD.Allgather(&number, 1, MPI::INT, vector.data(), 1, MPI::INT);
    else
        vector[0] = number;
}

void Parallel::allGather(double number, std::vector<double> &vector)
{
    if(isInitialized_)
        MPI::COMM_WORLD.Allgather(&number, 1, MPI::DOUBLE, vector.data(), 1, MPI::DOUBLE);
    else
        vector[0] = number;
}

void Parallel::barrier()
{
    MPI::COMM_WORLD.Barrier();
}

//*************************** Private Methods ************************************

void Parallel::loadBuffer(const std::vector<Vector3D> &vecs)
{
    if(commBuffer_.size() < 3*vecs.size())
        commBuffer_.resize(3*vecs.size());

    for(int i = 0, k = 0; i < vecs.size(); ++i, k += 3)
    {
        commBuffer_[k] = vecs[i].x;
        commBuffer_[k + 1] = vecs[i].y;
        commBuffer_[k + 2] = vecs[i].z;
    }
}

void Parallel::unloadBuffer(std::vector<Vector3D> &vecs)
{
    if(commBuffer_.size() < 3*vecs.size())
        commBuffer_.resize(3*vecs.size());

    for(int i = 0, k = 0; i < vecs.size(); ++i, k += 3)
    {
        vecs[i].x = commBuffer_[k];
        vecs[i].y = commBuffer_[k + 1];
        vecs[i].z = commBuffer_[k + 2];
    }
}

//***************************** Private Static Variables *******************************

std::vector<double> Parallel::commBuffer_;
int Parallel::procNo_ = 0, Parallel::nProcesses_ = 1;
bool Parallel::isInitialized_ = false;
