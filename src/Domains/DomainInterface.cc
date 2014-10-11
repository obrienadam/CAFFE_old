#include "DomainInterface.h"

DomainInterface::~DomainInterface()
{

}

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

void DomainInterface::allocate(Input &input)
{

    int nI = input.inputInts["nI"],
            nJ = input.inputInts["nJ"],
            nK = input.inputInts["nK"];

    for(scalarFieldMapItr_ = scalarFields_.begin(); scalarFieldMapItr_ != scalarFields_.end(); ++scalarFieldMapItr_)
    {

        scalarFieldMapItr_->second->allocate(nI, nJ, nK);

    }

    for(vectorFieldMapItr_ = vectorFields_.begin(); vectorFieldMapItr_ != vectorFields_.end(); ++vectorFieldMapItr_)
    {

        vectorFieldMapItr_->second->allocate(nI, nJ, nK);

    }

    for(tensorFieldMapItr_ = tensorFields_.begin(); tensorFieldMapItr_ != tensorFields_.end(); ++tensorFieldMapItr_)
    {

        tensorFieldMapItr_->second->allocate(nI, nJ, nK);

    }

}

DomainInterface& DomainInterface::computeTimeDerivative(SchemeInterface *scheme)
{



}

DomainInterface& DomainInterface::operator=(const DomainInterface& rhs)
{



}

DomainInterface& DomainInterface::operator+=(const DomainInterface& rhs)
{



}


DomainInterface& DomainInterface::operator*(const double rhs)
{



}
