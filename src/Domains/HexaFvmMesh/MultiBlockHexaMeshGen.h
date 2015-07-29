#ifndef MULTI_BLOCK_HEXA_MESH_GEN_H
#define MULTI_BLOCK_HEXA_MESH_GEN_H

<<<<<<< HEAD
#include "StructuredMesh.h"

class MultiBlockHexaMeshGen
{
public:

    MultiBlockHexaMeshGen();

    void generateMesh();
    void writeMeshFiles();

private:

    void getOwnershipRange(int nEntities, int nProcs, int procNo, int &iLower, int &iUpper, int &nEntitiesThisProc);

    StructuredMesh mesh_;
    int nBlocksI_, nBlocksJ_, nBlocksK_;
=======
class MultiBlockHexaMeshGen
{

>>>>>>> dd610b5763d66236723e62dc145c1660439db2a6
};

#endif
