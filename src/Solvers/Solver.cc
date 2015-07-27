#include "Solver.h"
#include "Output.h"

Solver::Solver(const Input &input, const HexaFvmMesh &mesh)
    :
      mesh_(mesh)
{
    indexMap_.initialize(mesh_.nCellsI(), mesh_.nCellsJ(), mesh_.nCellsK());

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
    int i;

    A_.resize(nMatrices);
    x_.resize(nVectors);
    b_.resize(nVectors);
    res_.resize(nVectors);

    for(i = 0; i < nMatrices; ++i)
    {
        A_[i].allocate(indexMap_.nActive(), indexMap_.nActive(), nnz);
    }

    for(i = 0; i < nVectors; ++i)
    {
        x_[i].allocate(indexMap_.nActive());
        b_[i].allocate(indexMap_.nActive());
        res_[i].allocate(indexMap_.nActive());
    }
}

void Solver::destroyMatrices()
{
    int i;

    for(i = 0; i < A_.size(); ++i)
    {
        A_[i].deallocate();
    }

    for(i = 0; i < x_.size(); ++i)
    {
        x_[i].deallocate();
        b_[i].deallocate();
        res_[i].deallocate();
    }
}

void Solver::zeroMatrices()
{
    int i;

    for(i = 0; i < A_.size(); ++i)
    {
        A_[i].zeroEntries();
    }

    for(i = 0; i < x_.size(); ++i)
    {
        x_[i].zeroEntries();
        b_[i].zeroEntries();
    }
}
