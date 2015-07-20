#include "Solver.h"
#include "Output.h"

Solver::Solver(const Input &input, const HexaFvmMesh &mesh)
    :
      mesh_(mesh),
      dE_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      dW_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      dN_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      dS_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      dT_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      dB_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      cE_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      cW_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      cN_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      cS_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      cT_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK()),
      cB_(mesh.nCellsI(), mesh.nCellsJ(), mesh.nCellsK())
{
    indexMap_.initialize(mesh_.nCellsI(), mesh_.nCellsJ(), mesh_.nCellsK());
    computeMeshMetrics();

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

void Solver::computeMeshMetrics()
{
    int i, j, k;

    for(k = 0; k < mesh_.nCellsK(); ++k)
    {
        for(j = 0; j < mesh_.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh_.nCellsI(); ++i)
            {
                dE_(i, j, k) = dot(mesh_.fAreaNormE(i, j, k), mesh_.fAreaNormE(i, j, k))/dot(mesh_.fAreaNormE(i, j, k), mesh_.rCellE(i, j, k));
                dW_(i, j, k) = dot(mesh_.fAreaNormW(i, j, k), mesh_.fAreaNormW(i, j, k))/dot(mesh_.fAreaNormW(i, j, k), mesh_.rCellW(i, j, k));
                dN_(i, j, k) = dot(mesh_.fAreaNormN(i, j, k), mesh_.fAreaNormN(i, j, k))/dot(mesh_.fAreaNormN(i, j, k), mesh_.rCellN(i, j, k));
                dS_(i, j, k) = dot(mesh_.fAreaNormS(i, j, k), mesh_.fAreaNormS(i, j, k))/dot(mesh_.fAreaNormS(i, j, k), mesh_.rCellS(i, j, k));
                dT_(i, j, k) = dot(mesh_.fAreaNormT(i, j, k), mesh_.fAreaNormT(i, j, k))/dot(mesh_.fAreaNormT(i, j, k), mesh_.rCellT(i, j, k));
                dB_(i, j, k) = dot(mesh_.fAreaNormB(i, j, k), mesh_.fAreaNormB(i, j, k))/dot(mesh_.fAreaNormB(i, j, k), mesh_.rCellB(i, j, k));

                cE_(i, j, k) = mesh_.fAreaNormE(i, j, k) - mesh_.rCellE(i, j, k)*dot(mesh_.fAreaNormE(i, j, k), mesh_.fAreaNormE(i, j, k))/dot(mesh_.fAreaNormE(i, j, k), mesh_.rCellE(i, j, k));
                cW_(i, j, k) = mesh_.fAreaNormW(i, j, k) - mesh_.rCellW(i, j, k)*dot(mesh_.fAreaNormW(i, j, k), mesh_.fAreaNormW(i, j, k))/dot(mesh_.fAreaNormW(i, j, k), mesh_.rCellW(i, j, k));
                cN_(i, j, k) = mesh_.fAreaNormN(i, j, k) - mesh_.rCellN(i, j, k)*dot(mesh_.fAreaNormN(i, j, k), mesh_.fAreaNormN(i, j, k))/dot(mesh_.fAreaNormN(i, j, k), mesh_.rCellN(i, j, k));
                cS_(i, j, k) = mesh_.fAreaNormS(i, j, k) - mesh_.rCellS(i, j, k)*dot(mesh_.fAreaNormS(i, j, k), mesh_.fAreaNormS(i, j, k))/dot(mesh_.fAreaNormS(i, j, k), mesh_.rCellS(i, j, k));
                cT_(i, j, k) = mesh_.fAreaNormT(i, j, k) - mesh_.rCellT(i, j, k)*dot(mesh_.fAreaNormT(i, j, k), mesh_.fAreaNormT(i, j, k))/dot(mesh_.fAreaNormT(i, j, k), mesh_.rCellT(i, j, k));
                cB_(i, j, k) = mesh_.fAreaNormB(i, j, k) - mesh_.rCellB(i, j, k)*dot(mesh_.fAreaNormB(i, j, k), mesh_.fAreaNormB(i, j, k))/dot(mesh_.fAreaNormB(i, j, k), mesh_.rCellB(i, j, k));
            }
        }
    }
}
