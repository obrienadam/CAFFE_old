#ifndef MULTI_BLOCK_HEXA_MESH_GEN_H
#define MULTI_BLOCK_HEXA_MESH_GEN_H

#include "StructuredMesh.h"

class MultiBlockHexaMeshGen
{
public:

    MultiBlockHexaMeshGen();

    void readFile();
    void generateMesh();
    void writeMeshFiles();

private:

    void getOwnershipRange(int nEntities, int nProcs, int procNo, int &iLower, int &iUpper, int &nEntitiesThisProc);

    StructuredMesh mesh_;
    int nBlocksI_, nBlocksJ_, nBlocksK_;

};

#endif
