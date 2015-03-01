#ifndef FV_SCHEME_H
#define FV_SCHEME_H

#include <vector>
#include <string>

#include "Field.h"
#include "HexaFvmMesh.h"

enum SchemeType{EXPLICIT, IMPLICIT};

class FvScheme
{
protected:

    std::string conservedFieldName_;
    HexaFvmMesh* meshPtr_;
    SchemeType schemeType;

public:

    FvScheme();

    virtual void initialize(HexaFvmMesh& mesh, std::string conservedFieldName = "phi");

    virtual void discretize() = 0;
    virtual void integrate(double timeStep) = 0;
};

#endif
