#include "DomainInterface.h"
#include "Output.h"

DomainInterface::~DomainInterface()
{

}

void DomainInterface::addField(ScalarField& scalarField)
{

    scalarFields_.push_back(&scalarField);

}

void DomainInterface::addField(VectorField& vectorField)
{

    vectorFields_.push_back(&vectorField);

}

void DomainInterface::addField(TensorField& tensorField)
{

    tensorFields_.push_back(&tensorField);

}

void DomainInterface::addAuxField(ScalarField& scalarField)
{

    scalarAuxFields_.push_back(&scalarField);

}

void DomainInterface::addAuxField(VectorField& vectorField)
{

    vectorAuxFields_.push_back(&vectorField);

}

void DomainInterface::addAuxField(TensorField& tensorField)
{

    tensorAuxFields_.push_back(&tensorField);

}

void DomainInterface::allocateFields(Input & input)
{

    int nI = input.inputInts["nI"],
            nJ = input.inputInts["nJ"],
            nK = input.inputInts["nK"],
            fieldNo;

    //- Initialize the main fields

    for(fieldNo = 0; fieldNo < scalarFields_.size(); ++fieldNo)
    {

        scalarFields_[fieldNo]->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized scalar field \"" + scalarFields_[fieldNo]->fieldName + "\".");

    }

    for(fieldNo = 0; fieldNo < vectorFields_.size(); ++fieldNo)
    {

        vectorFields_[fieldNo]->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized vector field \"" + vectorFields_[fieldNo]->fieldName + "\".");

    }

    for(fieldNo = 0; fieldNo < tensorFields_.size(); ++fieldNo)
    {

        tensorFields_[fieldNo]->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized tensor field \"" + tensorFields_[fieldNo]->fieldName + "\".");

    }

    //- Initialize the auxillary fields

    for(fieldNo = 0; fieldNo < scalarAuxFields_.size(); ++fieldNo)
    {

        scalarAuxFields_[fieldNo]->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized auxillary scalar field \"" + scalarAuxFields_[fieldNo]->fieldName + "\".");

    }

    for(fieldNo = 0; fieldNo < vectorAuxFields_.size(); ++fieldNo)
    {

        vectorAuxFields_[fieldNo]->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized auxillary vector field \"" + vectorAuxFields_[fieldNo]->fieldName + "\".");

    }

    for(fieldNo = 0; fieldNo < tensorAuxFields_.size(); ++fieldNo)
    {

        tensorAuxFields_[fieldNo]->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized auxillary tensor field \"" + tensorAuxFields_[fieldNo]->fieldName + "\".");

    }

    Output::printLine();
}

void DomainInterface::allocate(Input &input)
{

    allocateFields(input);

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
