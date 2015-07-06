#include <mpi.h>

#include "Parallel.h"

int Parallel::nProcesses()
{
    return MPI::COMM_WORLD.Get_size();
}

int Parallel::processNo()
{
    return MPI::COMM_WORLD.Get_rank();
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

void Parallel::sendMatrix(int root, int dest, Matrix &matrix)
{
    // Note - there is no check here to make sure the matrices on each process are the right size in the interest of efficiency
    if(processNo() == root)
        MPI::COMM_WORLD.Send(matrix.data(), matrix.nElements(), MPI::DOUBLE, root, dest);
    else if(processNo() == dest)
        MPI::COMM_WORLD.Recv(matrix.data(), matrix.nElements(), MPI::DOUBLE, root, dest);
}
