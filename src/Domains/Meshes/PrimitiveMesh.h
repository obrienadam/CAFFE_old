#ifndef PRIMITIVE_MESH_H
#define PRIMITIVE_MESH_H

#include "Point3D.h"
#include "SmartPointer3D.h"

enum Patch{BOUNDARY, INTERIOR};
enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};

class PrimitiveMesh
{

 protected:

  SmartPointer3D<Point3D> nodes_;

  Patch facePatches_[6];

 public:

  PrimitiveMesh();
  PrimitiveMesh(int nI, int nJ, int nK);
  ~PrimitiveMesh();

  void allocate(int nI, int nJ, int nK);

  Point3D& operator()(int i, int j, int k);
};

#endif // PRIMITIVEMESH_H
