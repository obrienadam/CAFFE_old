#ifndef HEXA_FDM_MESH_H
#define HEXA_FDM_MESH_H

#include "SmartPointer3D.h"
#include "Point3D.h"

class HexaFdmMesh
{
 private:

  SmartPointer3D<Point3D> nodes_;

 public:
  
  int eastPatch, westPatch, northPatch, southPatch, topPatch, bottomPatch;

  HexaFdmMesh();
  HexaFdmMesh(int nI, int nJ, int nK);

  Point3D& operator()(int i, int j, int k);
};

#endif
