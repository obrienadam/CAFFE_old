#ifndef HEXA_FDM_MESH_H
#define HEXA_FDM_MESH_H

#include "DomainInterface.h"
#include "StructuredMesh.h"

template <class STATE_TYPE>
class HexaFdmMesh : public StructuredMesh<STATE_TYPE>
{
 private:

 public:

  HexaFdmMesh();
  HexaFdmMesh(int nI, int nJ, int nK);

  void initialize(Input &input);

};

#include "HexaFdmMeshI.h"

#endif
