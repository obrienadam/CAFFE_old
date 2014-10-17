#ifndef HEXA_FDM_MESH_H
#define HEXA_FDM_MESH_H

#include "DomainInterface.h"
#include "StructuredMesh.h"

class HexaFdmMesh : public StructuredMesh
{
 private:

 public:

  HexaFdmMesh();
  HexaFdmMesh(int nI, int nJ, int nK);

  void allocate(Input &input);

};

#endif
