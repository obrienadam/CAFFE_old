#ifndef DOMAIN_INTERFACE_H
#define DOMAIN_INTERFACE_H

#include <map>
#include <cstdlib>

#include "ScalarField.h"
#include "VectorField.h"
#include "TensorField.h"


class DomainInterface
{

    typedef std::map<std::string, ScalarField*> ScalarFieldMap;
    typedef std::map<std::string, VectorField*> VectorFieldMap;
    typedef std::map<std::string, TensorField*> TensorFieldMap;


protected:

    ScalarFieldMap scalarFields_;
    VectorFieldMap vectorFields_;
    TensorFieldMap tensorFields_;

public:

    void addField(ScalarField& scalarField);
    void addField(VectorField& vectorField);
    void addField(TensorField& tensorField);

    ScalarField& scalar(const std::string& scalarFieldName);
    VectorField& vector(const std::string& vectorFieldName);
    TensorField& tensor(const std::string& tensorFieldName);
};

#endif
