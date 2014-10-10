#ifndef DOMAIN_INTERFACE_H
#define DOMAIN_INTERFACE_H

#include <map>
#include <cstdlib>

#include "ScalarField.h"
#include "VectorField.h"
#include "TensorField.h"
#include "Input.h"


class DomainInterface
{

protected:

    typedef std::map<std::string, ScalarField*> ScalarFieldMap;
    typedef std::map<std::string, VectorField*> VectorFieldMap;
    typedef std::map<std::string, TensorField*> TensorFieldMap;

    ScalarFieldMap scalarFields_;
    VectorFieldMap vectorFields_;
    TensorFieldMap tensorFields_;

public:

    virtual ~DomainInterface();

    void addField(ScalarField& scalarField);
    void addField(VectorField& vectorField);
    void addField(TensorField& tensorField);

    ScalarField& scalar(const std::string& scalarFieldName);
    VectorField& vector(const std::string& vectorFieldName);
    TensorField& tensor(const std::string& tensorFieldName);

    virtual void allocate(Input& input);
};

#endif
