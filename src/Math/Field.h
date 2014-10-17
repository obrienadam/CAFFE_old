#ifndef FIELD_H
#define FIELD_H

#include <string>

#include "Array3D.h"

template <class T>
class Field : public Array3D<T>
{

public:

    Field();

    Field(std::string fieldName, int nI = 0,
          int nJ = 0,
          int nK = 0);

    std::string fieldName;

    Field& operator+=(const Field& rhs);
    Field& operator-=(const Field& rhs);
    Field& operator*=(double rhs);
    Field& operator/=(double rhs);

    typedef typename Array3D<T>::iterator iterator;

};

#include "FieldI.h"

#endif
