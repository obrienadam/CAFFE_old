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
    virtual int nConservedVariables() = 0;

    virtual void discretize(std::vector<double>& timeDerivatives_) = 0;
    virtual void updateSolution(std::vector<double>& timeDerivatives_) = 0;
};

#endif
