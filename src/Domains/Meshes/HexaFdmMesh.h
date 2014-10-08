#ifndef HEXA_FDM_MESH_H
#define HEXA_FDM_MESH_H

#include "DomainInterface.h"
#include "PrimitiveMesh.h"

class HexaFdmMesh
        :
        public DomainInterface,
        public PrimitiveMesh
{
 private:

 public:

  HexaFdmMesh();
  HexaFdmMesh(int nI, int nJ, int nK);

  Point3D& operator()(int i, int j, int k);
};

#endif
