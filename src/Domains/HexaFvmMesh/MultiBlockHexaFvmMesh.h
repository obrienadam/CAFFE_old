#ifndef MULTI_BLOCK_HEXA_FVM_MESH_H
#define MULTI_BLOCK_HEXA_FVM_MESH_H

#include <vector>
#include <string>

#include "Input.h"
#include "HexaFvmMesh.h"

class MultiBlockHexaFvmMesh
{
public:
<<<<<<< HEAD
=======

    MultiBlockHexaFvmMesh();


>>>>>>> dd610b5763d66236723e62dc145c1660439db2a6

private:

<<<<<<< HEAD
    int nBlocks() const { return blocks_.size(); }

private:

=======
>>>>>>> dd610b5763d66236723e62dc145c1660439db2a6
    std::vector<HexaFvmMesh> blocks_;

    int nProcesses_, lProcess_, uProcess_;

};

#endif
