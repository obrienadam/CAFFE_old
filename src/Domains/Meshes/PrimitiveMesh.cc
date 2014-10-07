#include <cstdlib>
#include <algorithm>

#include "PrimitiveMesh.h"

PrimitiveMesh::PrimitiveMesh()
{
  std::fill_n(facePatches_, 6, INTERIOR);
}

PrimitiveMesh::PrimitiveMesh(int nI, int nJ, int nK)
  :
  PrimitiveMesh()
{
  allocate(nI, nJ, nK);
}

PrimitiveMesh::~PrimitiveMesh()
{

}

void PrimitiveMesh::allocate(int nI, int nJ, int nK)
{
  nodes_.allocate(nI, nJ, nK);
}

Point3D& PrimitiveMesh::operator()(int i, int j, int k)
{
    return nodes_(i, j, k);
}
