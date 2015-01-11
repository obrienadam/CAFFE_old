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

public:

    FvScheme();

    virtual void initialize(HexaFvmMesh& mesh, std::string conservedFieldName = "phi");
    virtual double computeFaceFlux(int i, int j, int k, Face face) = 0;
    virtual void computeSemiDiscreteForm() = 0;

};

#endif
