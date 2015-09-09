#include "Solver.h"
#include "Output.h"

Solver::Solver(const Input &input, const HexaFvmMesh &mesh)
    :
      mesh_(mesh)
{
    if(input.caseParameters.get<std::string>("Solver.timeAccurate") == "ON")
        solutionType_ = UNSTEADY;
    else
        solutionType_ = STEADY;

    Output::print("Solver", "Time accurate: " + input.caseParameters.get<std::string>("Solver.timeAccurate"));
}

void Solver::displayUpdateMessage()
{
    Output::print("Solver", "completed iteration.");
    Output::print("Solver", "BiCGStab iterations: " + std::to_string(biCGStabIters_));
    Output::printLine();
}

void Solver::createMatrices(int nMatrices, int nVectors, int nnz)
{
    A_.resize(nMatrices);
    x_.resize(nVectors);
    b_.resize(nVectors);
    res_.resize(nVectors);

    for(int i = 0; i < nMatrices; ++i)
    {
        A_[i].allocate(mesh_.iMap.nActive(), mesh_.iMap.nActive(), nnz);
    }

    for(int i = 0; i < nVectors; ++i)
    {
        x_[i].allocate(mesh_.iMap.nActive());
        b_[i].allocate(mesh_.iMap.nActive());
        res_[i].allocate(mesh_.iMap.nActive());
    }
}

void Solver::createMatrices(int nVariables, int nMatrices, int nVectors, int nnz)
{
    A_.resize(nMatrices);
    x_.resize(nVectors);
    b_.resize(nVectors);
    res_.resize(nVectors);

    for(int i = 0; i < nMatrices; ++i)
    {
        A_[i].allocate(nVariables*mesh_.iMap.nActive(), nVariables*mesh_.iMap.nActive(), nnz);
    }

    for(int i = 0; i < nVectors; ++i)
    {
        x_[i].allocate(nVariables*mesh_.iMap.nActive());
        b_[i].allocate(nVariables*mesh_.iMap.nActive());
        res_[i].allocate(nVariables*mesh_.iMap.nActive());
    }
}

void Solver::destroyMatrices()
{
    for(int i = 0; i < A_.size(); ++i)
    {
        A_[i].deallocate();
    }

    for(int i = 0; i < x_.size(); ++i)
    {
        x_[i].deallocate();
        b_[i].deallocate();
        res_[i].deallocate();
    }
}

void Solver::zeroMatrices()
{
    for(int i = 0; i < A_.size(); ++i)
    {
        A_[i].zeroEntries();
    }

    for(int i = 0; i < x_.size(); ++i)
    {
        x_[i].zeroEntries();
        b_[i].zeroEntries();
    }
}
