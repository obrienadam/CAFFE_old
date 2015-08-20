#ifndef MULTI_BLOCK_HEXA_MESH_GEN_H
#define MULTI_BLOCK_HEXA_MESH_GEN_H

#include "StructuredMesh.h"

class MultiBlockHexaMeshGen
{
public:

    MultiBlockHexaMeshGen();

    void initializeCartesianMesh(double xLength, double yLength, double zLength, int nCellsI, int nCellsJ, int nCellsK);

    void readFile(const std::string &directory);
    void generateMesh();
    void writeMeshFiles(const std::string &directory);

    StructuredMesh mesh;

private:

    void getOwnershipRange(int nEntities, int nProcs, int procNo, int &iLower, int &iUpper, int &nEntitiesThisProc);

    int nBlocksI_, nBlocksJ_, nBlocksK_;

};

#endif
