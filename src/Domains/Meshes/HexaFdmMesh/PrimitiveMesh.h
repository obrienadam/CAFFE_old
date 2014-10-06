#ifndef PRIMITIVE_MESH_H
#define PRIMITIVE_MESH_H

#include "Vector3D.h"
#include "SmartPointer3D.h"

class PrimitiveMesh
{

protected:

    SmartPointer3D<Vector3D> nodes_;

public:

    PrimitiveMesh();
    PrimitiveMesh(int nI, int nJ, int nK);

    Vector3D& operator()(int i, int j, int k) const;
};

#endif // PRIMITIVEMESH_H
