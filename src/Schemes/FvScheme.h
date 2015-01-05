#ifndef FV_SCHEME_H
#define FV_SCHEME_H

#include <vector>
#include <string>

#include "Field.h"
#include "HexaFvmMesh.h"

class FvScheme
{

protected:

    std::string conservedFieldName_;
    HexaFvmMesh* meshPtr_;
    Field<double> conservedFieldDerivatives_;

public:

    FvScheme();

    virtual void initialize(HexaFvmMesh& mesh, std::string conservedFieldName = "phi");

};

#endif
