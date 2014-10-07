#include <cstdlib>
#include <algorithm>

#include "PrimitiveMesh.h"

PrimitiveMesh::PrimitiveMesh()
{
  std::fill_n(facePatches_, 6, INTERIOR);
  boundaryMeshes_ = new PrimitiveMesh[6];
}

PrimitiveMesh::PrimitiveMesh(int nI, int nJ, int nK)
  :
  PrimitiveMesh()
{
  allocate(nI, nJ, nK);
}

PrimitiveMesh::~PrimitiveMesh()
{
  delete[] boundaryMeshes_;
}

void PrimitiveMesh::allocate(int nI, int nJ, int nK)
{
  nodes_.allocate(nI, nJ, nK);
}

Point3D& PrimitiveMesh::operator()(int i, int j, int k)
{
    return nodes_(i, j, k);
}

void PrimitiveMesh::createBoundary(Face location, Patch type, int order)
{
  if( location == EAST || location == WEST )
    {
      boundaryMeshes_[location].allocate(order, nodes_.sizeJ(), nodes_.sizeK());
    }
  else if (location == NORTH || location == BOTTOM)
    {
      boundaryMeshes_[location].allocate(nodes_.sizeI(), order, nodes_.sizeK());
    }
  else if (location == TOP || location == BOTTOM)
    {
      boundaryMeshes_[location].allocate(nodes_.sizeI(), nodes_.sizeJ(), order);
    }
  else
    {
      throw "Invalid face location specified in PrimitiveMesh::createBoundary.";
    }

  facePatches_[location] = type;
}
