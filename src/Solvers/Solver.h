#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

#include "Input.h"
#include "HexaFvmMesh.h"
#include "SparseMatrix.h"
#include "SparseVector.h"
#include "IndexMap.h"
#include "Time.h"

class Solver
{
public:

    enum SolutionType{STEADY, UNSTEADY};

    Solver(const Input &input, const HexaFvmMesh &mesh);

    virtual double solve(double timeStep) = 0;
    virtual void displayUpdateMessage();

protected:

    void createMatrices(int nMatrices, int nVectors, int nnz);
    void destroyMatrices();
    void zeroMatrices();
    void computeMeshMetrics();

    const HexaFvmMesh &mesh_;

    SolutionType solutionType_;
    std::vector<SparseMatrix> A_;
    std::vector<SparseVector> x_, b_, res_;
    IndexMap indexMap_;

    Array3D<double> dE_, dW_, dN_, dS_, dT_, dB_;
    Array3D<Vector3D> cE_, cW_, cN_, cS_, cT_, cB_;
    Time time_;

    int biCGStabIters_;

    friend class InitialConditions;
};

#endif
