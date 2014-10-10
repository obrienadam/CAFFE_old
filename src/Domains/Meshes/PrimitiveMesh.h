#ifndef PRIMITIVE_MESH_H
#define PRIMITIVE_MESH_H

#include "DomainInterface.h"
#include "Point3D.h"
#include "SmartPointer3D.h"

enum Patch{BOUNDARY, INTERIOR};
enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};

class PrimitiveMesh : public DomainInterface
{

 protected:

  SmartPointer3D<Point3D> nodes_;

  Patch facePatches_[6];

 public:

  PrimitiveMesh();
  virtual ~PrimitiveMesh();

  virtual void allocate(Input& input);

};

#endif // PRIMITIVEMESH_H
