#ifndef DOMAIN_INTERFACE_H
#define DOMAIN_INTERFACE_H

#include <map>
#include <cstdlib>

#include "ScalarField.h"
#include "FieldInterface.h"

class DomainInterface
{

 protected:

 public:

    std::map<std::string, FieldInterface*> fields;

    void addField(FieldInterface* field)
    {

        fields[field->fieldName] = field;

    }

};

#endif
