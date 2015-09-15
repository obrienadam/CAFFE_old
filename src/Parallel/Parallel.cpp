#include "Parallel.h"

const int Parallel::PROC_NULL = MPI::PROC_NULL;

/********************** Public methods ***********************/

void Parallel::initialize()
{
    MPI::Init();
    procNo_ = MPI::COMM_WORLD.Get_rank();
    nProcesses_ = MPI::COMM_WORLD.Get_size();
    requests_.reserve(100);
    isInitialized_ = true;
}

void Parallel::finalize()
{
    MPI::Finalize();
    procNo_ = 0;
    nProcesses_ = 1;
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
        MPI::COMM_WORLD.Bcast(buffers_.initCollectiveBuffer(Parallel::processNo() == source, vecs), 3*vecs.size(), MPI::DOUBLE, source);
        buffers_.freeLastVector3DRecvBuffer();
    }
}

int Parallel::sum(int number)
{
    int sum = number;

    if(isInitialized_)
        MPI::COMM_WORLD.Allreduce(&number, &sum, 1, MPI::INT, MPI::SUM);

    return sum;
}

double Parallel::sum(double number)
{
    double sum = number;

    if(isInitialized_)
        MPI::COMM_WORLD.Allreduce(&number, &sum, 1, MPI::DOUBLE, MPI::SUM);

    return sum;
}

int Parallel::min(int number)
{
    int min = number;

    if(isInitialized_)
        MPI::COMM_WORLD.Allreduce(&number, &min, 1, MPI::INT, MPI::MIN);

    return min;
}

double Parallel::min(double number)
{
    double min = number;

    if(isInitialized_)
        MPI::COMM_WORLD.Allreduce(&number, &min, 1, MPI::DOUBLE, MPI::MIN);

    return min;
}

int Parallel::max(int number)
{
    int max = number;

    if(isInitialized_)
        MPI::COMM_WORLD.Allreduce(&number, &max, 1, MPI::INT, MPI::MAX);

    return max;
}

double Parallel::max(double number)
{
    double max = number;

    if(isInitialized_)
        MPI::COMM_WORLD.Allreduce(&number, &max, 1, MPI::DOUBLE, MPI::MAX);

    return max;
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
            MPI::COMM_WORLD.Send(buffers_.initSendBuffer(vecs), 3*vecs.size(), MPI::DOUBLE, dest, source);
            buffers_.freeLastVector3DSendBuffer();
        }
        if(processNo() == dest)
        {
            MPI::COMM_WORLD.Recv(buffers_.initRecvBuffer(vecs), 3*vecs.size(), MPI::DOUBLE, source, source);
            buffers_.freeLastVector3DRecvBuffer();
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

void Parallel::iSend(int source, int dest, int tag, const std::vector<int> &ints)
{
    if(isInitialized_)
    {
        if(processNo() == source)
            requests_.push_back(MPI::COMM_WORLD.Isend(ints.data(), ints.size(), MPI::INT, dest, tag));
    }
}

void Parallel::iRecv(int source, int dest, int tag, std::vector<int> &ints)
{
    if(isInitialized_)
    {
        if(processNo() == dest)
            requests_.push_back(MPI::COMM_WORLD.Irecv(ints.data(), ints.size(), MPI::INT, source, tag));
    }
}

void Parallel::iSend(int source, int dest, int tag, const std::vector<double> &doubles)
{
    if(isInitialized_)
    {
        if(processNo() == source)
            requests_.push_back(MPI::COMM_WORLD.Isend(doubles.data(), doubles.size(), MPI::DOUBLE, dest, tag));
    }
}

void Parallel::iRecv(int source, int dest, int tag, std::vector<double> &doubles)
{
    if(isInitialized_)
    {
        if(processNo() == dest)
            requests_.push_back(MPI::COMM_WORLD.Irecv(buffers_.initRecvBuffer(doubles), doubles.size(), MPI::DOUBLE, source, tag));
    }
}

void Parallel::iSend(int source, int dest, int tag, const std::vector<Vector3D> &vecs)
{
    if(isInitialized_ && processNo() == source)
        requests_.push_back(MPI::COMM_WORLD.Isend(buffers_.initSendBuffer(vecs), 3*vecs.size(), MPI::DOUBLE, dest, tag));
}

void Parallel::iRecv(int source, int dest, int tag, std::vector<Vector3D> &vecs)
{
    if(isInitialized_ && processNo() == dest)
        requests_.push_back(MPI::COMM_WORLD.Irecv(buffers_.initRecvBuffer(vecs), 3*vecs.size(), MPI::DOUBLE, source, tag));
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
        buffers_.freeAllBuffers();
    }
}

//***************************** Private Static Variables *******************************

ParallelBuffers Parallel::buffers_;
std::vector<MPI::Request> Parallel::requests_;
int Parallel::procNo_ = 0, Parallel::nProcesses_ = 1;
bool Parallel::isInitialized_ = false;
