#ifndef MULTI_BLOCK_HEXA_FVM_MESH_H
#define MULTI_BLOCK_HEXA_FVM_MESH_H

#include <vector>
#include "HexaFvmMesh.h"

class MultiBlockHexaFvmMesh
{

private:

    std::vector<HexaFvmMesh> blocks_;

public:

    MultiBlockHexaFvmMesh(int nBlocks = 0);

    void allocate(int nBlocks);
};

#endif
