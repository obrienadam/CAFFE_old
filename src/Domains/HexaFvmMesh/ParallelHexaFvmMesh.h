#ifndef PARALLEL_HEXA_FVM_MESH_H
#define PARALLEL_HEXA_FVM_MESH_H

#include <vector>
#include <string>
#include <memory>
#include <array>

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

    virtual std::string meshStats();

    int nSubDomains() const { return Parallel::nProcesses(); }

    std::shared_ptr< std::array<int, 6> > getAdjProcNoPtr() const { return adjProcNoPtr_; }

    void writeTec360(double time, const std::string &directory);
    void writeBoundaryMeshes(double time, const std::string &directory);

    virtual void changeName(const std::string &newName);

private:

    void initializeSubDomains(const StructuredMesh &tmpMesh, int nSubDomainsI, int nSubDomainsJ, int nSubDomainsK);

    std::shared_ptr< std::array<int, 6> > adjProcNoPtr_;
};

#endif
