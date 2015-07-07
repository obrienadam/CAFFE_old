#include "MultiBlockHexaFvmMesh.h"

MultiBlockHexaFvmMesh::MultiBlockHexaFvmMesh(int nBlocks)
{
    allocate(nBlocks);
}

void MultiBlockHexaFvmMesh::allocate(int nBlocks)
{
    blocks_.resize(nBlocks);
}
