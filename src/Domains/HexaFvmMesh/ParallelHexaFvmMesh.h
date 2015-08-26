#ifndef PARALLEL_HEXA_FVM_MESH_H
#define PARALLEL_HEXA_FVM_MESH_H

#include <vector>
#include <string>

#include "Input.h"
#include "HexaFvmMesh.h"
#include "Parallel.h"

class ParallelHexaFvmMesh : public HexaFvmMesh
{
public:

    ParallelHexaFvmMesh();

    virtual void initialize(const std::string &filename);
    virtual void initialize(const Array3D<Point3D> &nodes);
    virtual void initializeCartesianMesh(double xLength, double yLength, double zLength, int nCellsI, int nCellsJ, int nCellsK);

    int nSubDomains() const { return Parallel::nProcesses(); }

    const int* adjacentSubDomainProcNo() { return adjacentSubDomainProcNo_; }

    void writeTec360(double time, const std::string &directory);

private:

    void initializeSubDomains(const StructuredMesh &tmpMesh);

    std::vector<HexaFvmMesh> boundaryMeshes_;
    int adjacentSubDomainProcNo_[6];
};

#endif
