#ifndef PRIMITIVE_MESH_H
#define PRIMITIVE_MESH_H

#include "Point3D.h"
#include "SmartPointer3D.h"

class PrimitiveMesh
{

protected:

    SmartPointer3D<Point3D> nodes_;

public:

    PrimitiveMesh();
    PrimitiveMesh(int nI, int nJ, int nK);

    Vector3D& operator()(int i, int j, int k);
};

#endif // PRIMITIVEMESH_H
