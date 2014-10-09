#ifndef FIELD_H
#define FIELD_H

#include <string>

#include "SmartPointer3D.h"

template <class T>
class Field : public SmartPointer3D<T>
{

public:

    Field();

    Field(std::string fieldName, int nI, int nJ, int nK);

    std::string fieldName;

    Field& operator+=(const Field& rhs);
    Field& operator-=(const Field& rhs);
    Field& operator*=(double rhs);
    Field& operator/=(double rhs);

};

#include "FieldI.h"

#endif
