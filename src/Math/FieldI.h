#ifndef FIELD_I_H
#define FIELD_I_H

#include "Field.h"

template <class T>
Field<T>::Field()
{

}

template <class T>
Field<T>::Field(std::string fieldName,
                int nI,
                int nJ,
                int nK)
    :
      fieldName(fieldName),
      Array3D<T>(nI, nJ, nK)
{

}

#endif
