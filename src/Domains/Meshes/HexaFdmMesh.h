#ifndef HEXA_FDM_MESH_H
#define HEXA_FDM_MESH_H

#include "DomainInterface.h"
#include "PrimitiveMesh.h"

class HexaFdmMesh : public PrimitiveMesh
{
 private:

 public:

  HexaFdmMesh();
  HexaFdmMesh(int nI, int nJ, int nK);

};

#endif
