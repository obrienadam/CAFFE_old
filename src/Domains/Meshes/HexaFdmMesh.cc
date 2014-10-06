#include "HexaFdmMesh.h"

HexaFdmMesh::HexaFdmMesh()
{

}

HexaFdmMesh::HexaFdmMesh(int nI, int nJ, int nK)
:
  PrimitiveMesh(nI, nJ, nK)
{
  
}

Point3D& HexaFdmMesh::operator ()(int i, int j, int k)
{
    return nodes_(i, j, k);
}
