#ifndef FV_SCHEME_H
#define FV_SCHEME_H

#include <vector>
#include <string>

#include "HexaFvmMesh.h"

class FvScheme
{

protected:

    HexaFvmMesh* mesh_;

public:

    FvScheme();

    virtual void setMeshPointer(HexaFvmMesh* mesh);

};

#endif
