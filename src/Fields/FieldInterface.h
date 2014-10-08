#ifndef FIELD_INTERFACE_H
#define FIELD_INTERFACE_H

#include <string>

class FieldInterface
{

public:

    FieldInterface(){}

    FieldInterface(std::string fieldName)
    :
      fieldName(fieldName)
    {

    }

    std::string fieldName;

};

#endif // FIELDINTERFACE_H
