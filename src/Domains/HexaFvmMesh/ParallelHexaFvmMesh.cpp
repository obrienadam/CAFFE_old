#include "ParallelHexaFvmMesh.h"

ParallelHexaFvmMesh::ParallelHexaFvmMesh()
    :
      boundaryMeshes_(6)
{

}

void ParallelHexaFvmMesh::initialize(const std::string &filename)
{
    StructuredMesh tmpMesh;

    if(Parallel::isMainProcessor())
    {
        tmpMesh.initialize(filename);
    }

    initializeSubDomains(tmpMesh);
}

void ParallelHexaFvmMesh::initialize(const Array3D<Point3D> &nodes)
{
    StructuredMesh tmpMesh;

    if(Parallel::isMainProcessor())
    {
        tmpMesh.initialize(nodes);
    }

    initializeSubDomains(tmpMesh);
}

void ParallelHexaFvmMesh::initializeCartesianMesh(double xLength, double yLength, double zLength, int nCellsI, int nCellsJ, int nCellsK)
{
    StructuredMesh tmpMesh;

    if(Parallel::isMainProcessor())
    {
        tmpMesh.initializeCartesianMesh(xLength, yLength, zLength, nCellsI, nCellsJ, nCellsK);
    }

    initializeSubDomains(tmpMesh);
}

void ParallelHexaFvmMesh::writeTec360(double time, const std::string &directoryName)
{

}

//************************* Private methods **********************************

void ParallelHexaFvmMesh::initializeSubDomains(const StructuredMesh &tmpMesh)
{
    int nCellsIGlobal, nCellsJGlobal, nCellsKGlobal;
    std::vector<int> nCellsILocal, nCellsJLocal, nCellsKLocal;
    int nSubDomainsI, nSubDomainsJ, nSubDomainsK, subDomainNo;
    std::vector<Point3D> vertices(8);

    Output::issueWarning("ParallelHexaFvmMesh", "initializeSubDomains", "number of sub-domains currently fixed. May not be suitable for all computations.");
    nSubDomainsI = 2;
    nSubDomainsJ = 2;
    nSubDomainsK = 2;

    nCellsIGlobal = Parallel::broadcast(tmpMesh.nNodesI() + 1, Parallel::mainProcNo());
    nCellsJGlobal = Parallel::broadcast(tmpMesh.nNodesJ() + 1, Parallel::mainProcNo());
    nCellsKGlobal = Parallel::broadcast(tmpMesh.nNodesK() + 1, Parallel::mainProcNo());

    nCellsILocal.resize(Parallel::nProcesses());
    nCellsJLocal.resize(Parallel::nProcesses());
    nCellsKLocal.resize(Parallel::nProcesses());

    Parallel::allGather(nCellsIGlobal/nSubDomainsI, nCellsILocal);
    Parallel::allGather(nCellsJGlobal/nSubDomainsJ, nCellsJLocal);
    Parallel::allGather(nCellsKGlobal/nSubDomainsK, nCellsKLocal);

    if(Parallel::isMainProcessor())
    {
        vertices[0] = tmpMesh.node(0, 0, 0);
        vertices[1] = tmpMesh.node(tmpMesh.uNodeI(), 0, 0);
        vertices[2] = tmpMesh.node(tmpMesh.uNodeI(), tmpMesh.uNodeJ(), 0);
        vertices[3] = tmpMesh.node(0, tmpMesh.uNodeJ(), 0);
        vertices[4] = tmpMesh.node(0, 0, tmpMesh.uNodeK());
        vertices[5] = tmpMesh.node(tmpMesh.uNodeI(), 0, tmpMesh.uNodeK());
        vertices[6] = tmpMesh.node(tmpMesh.uNodeI(), tmpMesh.uNodeJ(), tmpMesh.uNodeK());
        vertices[7] = tmpMesh.node(0, tmpMesh.uNodeJ(), tmpMesh.uNodeK());

        HexaFvmMesh::initialize(vertices, nSubDomainsI, nSubDomainsJ, nSubDomainsK);
    }

    for(int k = 0; k < nSubDomainsK; ++k)
    {
        for(int j = 0; j < nSubDomainsJ; ++j)
        {
            for(int i = 0; i < nSubDomainsI; ++i)
            {
                subDomainNo = k*nSubDomainsI*nSubDomainsJ + j*nSubDomainsI + i;

                if(Parallel::isMainProcessor())
                {
                    vertices[0] = tmpMesh.node(i, j, k);
                    vertices[1] = tmpMesh.node(i + 1, j, k);
                    vertices[2] = tmpMesh.node(i + 1, j + 1, k);
                    vertices[3] = tmpMesh.node(i, j + 1, k);
                    vertices[4] = tmpMesh.node(i, j, k + 1);
                    vertices[5] = tmpMesh.node(i + 1, j, k + 1);
                    vertices[6] = tmpMesh.node(i + 1, j + 1, k + 1);
                    vertices[7] = tmpMesh.node(i, j + 1, k + 1);
                }

                Parallel::send(Parallel::mainProcNo(), subDomainNo, vertices);
            }
        }
    }
}
