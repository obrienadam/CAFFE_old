#include "DomainInterface.h"

void DomainInterface::addField(ScalarField& scalarField)
{

    scalarFields_[scalarField.fieldName] = &scalarField;

}

void DomainInterface::addField(VectorField& vectorField)
{

    vectorFields_[vectorField.fieldName] = &vectorField;

}

void DomainInterface::addField(TensorField& tensorField)
{

    tensorFields_[tensorField.fieldName] = &tensorField;

}

ScalarField& DomainInterface::scalar(const std::string &scalarFieldName)
{

    return *(scalarFields_[scalarFieldName]);

}

VectorField& DomainInterface::vector(const std::string &vectorFieldName)
{

    return *(vectorFields_[vectorFieldName]);

}

TensorField& DomainInterface::tensor(const std::string &tensorFieldName)
{

    return *(tensorFields_[tensorFieldName]);

}
