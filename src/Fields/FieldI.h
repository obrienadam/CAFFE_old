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

template <class T>
Field<T>& Field<T>::operator+=(const Field<T>& rhs)
{



}

template <class T>
Field<T>& Field<T>::operator-=(const Field<T>& rhs)
{



}

template <class T>
Field<T>& Field<T>::operator*=(double rhs)
{



}

template <class T>
Field<T>& Field<T>::operator/=(double rhs)
{



}

#endif
