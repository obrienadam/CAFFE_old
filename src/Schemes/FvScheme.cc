#include <cstdlib>

#include "FvScheme.h"
#include "Output.h"

FvScheme::FvScheme()
    :
      conservedFieldName_("phi"),
      meshPtr_(NULL)
{

}

void FvScheme::initialize(HexaFvmMesh &mesh, std::string conservedFieldName)
{
    meshPtr_ = &mesh;
    conservedFieldName_ = conservedFieldName;
}
