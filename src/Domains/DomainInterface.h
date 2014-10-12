#ifndef DOMAIN_INTERFACE_H
#define DOMAIN_INTERFACE_H

#include <vector>

#include "ScalarField.h"
#include "VectorField.h"
#include "TensorField.h"
#include "Input.h"
#include "SchemeInterface.h"


class DomainInterface
{

protected:

    typedef std::vector<ScalarField*> ScalarFieldVec;
    typedef std::vector<VectorField*> VectorFieldVec;
    typedef std::vector<TensorField*> TensorFieldVec;

    ScalarFieldVec scalarFields_;
    VectorFieldVec vectorFields_;
    TensorFieldVec tensorFields_;

    ScalarFieldVec scalarAuxFields_;
    VectorFieldVec vectorAuxFields_;
    TensorFieldVec tensorAuxFields_;

public:

    virtual ~DomainInterface();

    void addField(ScalarField& scalarField);
    void addField(VectorField& vectorField);
    void addField(TensorField& tensorField);

    void addAuxField(ScalarField& scalarField);
    void addAuxField(VectorField& vectorField);
    void addAuxField(TensorField& tensorField);

    void allocateFields(Input& input);
    virtual void allocate(Input& input);

};

#endif
