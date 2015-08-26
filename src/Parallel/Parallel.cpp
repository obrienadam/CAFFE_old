#include "Parallel.h"

/********************** Public methods ***********************/

void Parallel::initialize()
{
    MPI::Init();
    procNo_ = MPI::COMM_WORLD.Get_rank();
    nProcesses_ = MPI::COMM_WORLD.Get_size();
    commBuffer_.resize(100000);
    requests_.reserve(10000);
    isInitialized_ = true;
}

void Parallel::finalize()
{
    MPI::Finalize();
    procNo_ = 0;
    nProcesses_ = 1;
    commBuffer_.erase(commBuffer_.begin(), commBuffer_.end());
    requests_.erase(requests_.begin(), requests_.end());
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
    if(isInitialized_)
        MPI::COMM_WORLD.Bcast(&number, 1, MPI::INT, source);

    return number;
}

double Parallel::broadcast(double number, int source)
{
    if(isInitialized_)
        MPI::COMM_WORLD.Bcast(&number, 1, MPI::DOUBLE, source);

    return number;
}

Vector3D Parallel::broadcast(const Vector3D &vec, int source)
{
    double buff[3] = {vec.x, vec.y, vec.z};

    if(isInitialized_)
        MPI::COMM_WORLD.Bcast(buff, 3, MPI::DOUBLE, source);

    return Vector3D(buff[0], buff[1], buff[2]);
}

void Parallel::broadcast(int source, std::vector<double> &doubles)
{
    if(isInitialized_)
        MPI::COMM_WORLD.Bcast(doubles.data(), doubles.size(), MPI::DOUBLE, source);
}

void Parallel::broadcast(int source, std::vector<Vector3D> &vecs)
{
    if(isInitialized_)
    {
        if(processNo() == source)
            loadBuffer(vecs);

        MPI::COMM_WORLD.Bcast(commBuffer_.data(), 3*vecs.size(), MPI::DOUBLE, source);

        if(processNo() != source)
            unloadBuffer(vecs);
    }
}

void Parallel::send(int source, int dest, std::vector<double> &doubles)
{
    if(isInitialized_)
    {
        if(processNo() == source)
            MPI::COMM_WORLD.Send(doubles.data(), doubles.size(), MPI::DOUBLE, dest, source);
        if(processNo() == dest)
            MPI::COMM_WORLD.Recv(doubles.data(), doubles.size(), MPI::DOUBLE, source, source);
    }
}

void Parallel::send(int source, int dest, std::vector<Vector3D> &vecs)
{
    if(isInitialized_)
    {
        if(processNo() == source)
        {
            loadBuffer(vecs);
            MPI::COMM_WORLD.Send(commBuffer_.data(), 3*vecs.size(), MPI::DOUBLE, dest, source);
        }
        if(processNo() == dest)
        {
            MPI::COMM_WORLD.Recv(commBuffer_.data(), 3*vecs.size(), MPI::DOUBLE, source, source);
            unloadBuffer(vecs);
        }
    }
}

void Parallel::send(int source, int dest, const Array3D<int> &sourceArray3D, Array3D<int> &destArray3D)
{
    int dimensions[3] = {sourceArray3D.sizeI(), sourceArray3D.sizeJ(), sourceArray3D.sizeK()};

    if(isInitialized_)
    {
        if(processNo() == source)
        {
            MPI::COMM_WORLD.Send(dimensions, 3, MPI::INT, dest, source);
            MPI::COMM_WORLD.Send(sourceArray3D.data(), sourceArray3D.size(), MPI::INT, dest, source);
        }
        if(processNo() == dest)
        {
            MPI::COMM_WORLD.Recv(dimensions, 3, MPI::INT, source, source);
            destArray3D.resize(dimensions[0], dimensions[1], dimensions[2]);
            MPI::COMM_WORLD.Recv(destArray3D.data(), destArray3D.size(), MPI::INT, source, source);
        }
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
    if(isInitialized_)
        MPI::COMM_WORLD.Barrier();
}

void Parallel::waitAll()
{
    if(isInitialized_)
    {
        MPI::Request::Waitall(requests_.size(), requests_.data());
        requests_.clear();
    }
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
std::vector<MPI::Request> Parallel::requests_;
int Parallel::procNo_ = 0, Parallel::nProcesses_ = 1;
bool Parallel::isInitialized_ = false;
