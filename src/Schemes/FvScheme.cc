#include "FvScheme.h"
#include "Output.h"

FvScheme::FvScheme()
    :
      conservedFieldName_("phi")
{

}

void FvScheme::initialize(HexaFvmMesh &mesh, std::string conservedFieldName)
{

    meshPtr_ = &mesh;
    conservedFieldName_ = conservedFieldName;

}
