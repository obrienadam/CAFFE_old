#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include "DomainInterface.h"
#include "Point3D.h"
#include "Array3D.h"

enum Patch{BOUNDARY, INTERIOR};
enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};

class StructuredMesh : public DomainInterface
{

 protected:

  Array3D<Point3D> nodes_;

  Patch facePatches_[6];

 public:

  StructuredMesh();
  ~StructuredMesh();

  virtual void allocate(Input& input);
  virtual int size();

};

#endif // PRIMITIVEMESH_H
