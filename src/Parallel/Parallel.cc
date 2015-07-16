#include <mpi.h>

#include "Parallel.h"

int Parallel::commBufferSize_;
std::vector<double> Parallel::commBuffer_;

/********************** Public methods ***********************/

int Parallel::nProcesses()
{
    return MPI::COMM_WORLD.Get_size();
}

int Parallel::processNo()
{
    return MPI::COMM_WORLD.Get_rank();
}

bool Parallel::isThisProcessor(int source)
{
    if(source == processNo())
        return true;
    else
        return false;
}

bool Parallel::isMainProcessor()
{
    if(0 == processNo())
        return true;
    else
        return false;
}

void Parallel::ownershipRange(int nEntities, int &iLower, int &iUpper, int &nEntitiesThisProc)
{
    int nRemainingEntities;

    nEntitiesThisProc = nEntities/nProcesses();
    nRemainingEntities = nEntities%nProcesses();

    if(processNo() < nRemainingEntities)
    {
        ++nEntitiesThisProc;
        iLower = processNo()*nEntitiesThisProc;
        iUpper = iLower + nEntitiesThisProc - 1;
    }
    else
    {
        iLower = (nEntitiesThisProc + 1)*nRemainingEntities + (processNo() - nRemainingEntities)*nEntitiesThisProc;
        iUpper = iLower + nEntitiesThisProc - 1;
    }
}

int Parallel::sum(int number)
{
    MPI::COMM_WORLD.Allreduce(&number, &number, 1, MPI::INT, MPI::SUM);

    return number;
}

int Parallel::min(int number)
{
    MPI::COMM_WORLD.Allreduce(&number, &number, 1, MPI::INT, MPI::MIN);

    return number;
}

int Parallel::max(int number)
{
    MPI::COMM_WORLD.Allreduce(&number, &number, 1, MPI::INT, MPI::MAX);

    return number;
}

double Parallel::sum(double number)
{
    MPI::COMM_WORLD.Allreduce(&number, &number, 1, MPI::DOUBLE, MPI::SUM);

    return number;
}

double Parallel::min(double number)
{
    MPI::COMM_WORLD.Allreduce(&number, &number, 1, MPI::DOUBLE, MPI::MIN);

    return number;
}

double Parallel::max(double number)
{
    MPI::COMM_WORLD.Allreduce(&number, &number, 1, MPI::DOUBLE, MPI::MAX);

    return number;
}

void Parallel::send(int source, int dest, std::vector<double> &vector)
{
    if(processNo() == source)
        MPI::COMM_WORLD.Send(vector.data(), vector.size(), MPI::DOUBLE, dest, 0);
    if(processNo() == dest)
        MPI::COMM_WORLD.Recv(vector.data(), vector.size(), MPI::DOUBLE, source, 0);
}

void Parallel::send(int source, int dest, Array3D<double> &doubleArray3D)
{
    if(processNo() == source)
        MPI::COMM_WORLD.Send(doubleArray3D.data(), doubleArray3D.size(), MPI::DOUBLE, dest, 0);
    if(processNo() == dest)
        MPI::COMM_WORLD.Recv(doubleArray3D.data(), doubleArray3D.size(), MPI::DOUBLE, source, 0);
}

void Parallel::send(int source, int dest, Array3D<Vector3D> &vector3DArray3D)
{
    int i, k, size = 3*vector3DArray3D.size();

    if(commBufferSize_ < size)
    {
        commBufferSize_ = size;
        commBuffer_.resize(commBufferSize_);
    }

    if(processNo() == source)
    {
        for(i = 0; i < size; i += 3)
        {
            k = i/3;
            commBuffer_[i] = vector3DArray3D(k).x;
            commBuffer_[i + 1] = vector3DArray3D(k).y;
            commBuffer_[i + 2] = vector3DArray3D(k).z;
        }

        MPI::COMM_WORLD.Send(commBuffer_.data(), size, MPI::DOUBLE, dest, 0);
    }
    if(processNo() == dest)
    {
        MPI::COMM_WORLD.Recv(commBuffer_.data(), size, MPI::DOUBLE, source, 0);

        for(i = 0; i < size; i += 3)
        {
            k=i/3;
            vector3DArray3D(k).x = commBuffer_[i];
            vector3DArray3D(k).y = commBuffer_[i + 1];
            vector3DArray3D(k).z = commBuffer_[i + 2];
        }
    }
}

void Parallel::send(int source, int dest, Matrix &matrix)
{
    // Note - there is no check here to make sure the matrices on each process are the right size in the interest of efficiency
    if(processNo() == source)
        MPI::COMM_WORLD.Send(matrix.data(), matrix.nElements(), MPI::DOUBLE, dest, 0);
    if(processNo() == dest)
        MPI::COMM_WORLD.Recv(matrix.data(), matrix.nElements(), MPI::DOUBLE, source, 0);
}

void Parallel::barrier()
{
    MPI::COMM_WORLD.Barrier();
}

/********************** Private methods ***********************/

void Parallel::initialize()
{
    MPI::Init();

    commBufferSize_ = 1000000;
    commBuffer_.resize(commBufferSize_);
}

void Parallel::finalize()
{
    MPI::Finalize();
}
