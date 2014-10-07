#ifndef HEXA_FDM_MESH_H
#define HEXA_FDM_MESH_H

#include "PrimitiveMesh.h"

class HexaFdmMesh : public PrimitiveMesh
{
 private:

 public:

  HexaFdmMesh();
  HexaFdmMesh(int nI, int nJ, int nK);

  Point3D& operator()(int i, int j, int k);
};

#endif
