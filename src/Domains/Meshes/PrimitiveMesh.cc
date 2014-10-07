#include "PrimitiveMesh.h"

PrimitiveMesh::PrimitiveMesh()
{

}

PrimitiveMesh::PrimitiveMesh(int nI, int nJ, int nK)
    :
      nodes_(nI, nJ, nK)
{

}

Vector3D& PrimitiveMesh::operator()(int i, int j, int k)
{
    return nodes_(i, j, k);
}
