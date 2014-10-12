#ifndef SMART_POINTER_I_H
#define SMART_POINTER_I_H

#include <cstdlib>

#include "SmartPointer.h"

template <class T>
SmartPointer<T>::SmartPointer()
    :
      n_(0),
      data_(NULL),
      isOriginalObject_(true)
{



}

template <class T>
SmartPointer<T>::SmartPointer(int n)
    :
      SmartPointer<T>()
{

    allocate(n);

}

template <class T>
SmartPointer<T>::SmartPointer(T *ptr)
    :
      n_(1),
      data_(ptr),
      isOriginalObject_(false)
{

}

template <class T>
SmartPointer<T>::SmartPointer(const SmartPointer<T>& other)
    :
      n_(other.n_),
      data_(other.data_),
      isOriginalObject_(false)
{

}

template <class T>
SmartPointer<T>::~SmartPointer()
{

    if(isOriginalObject_)
        deallocate();

}

template <class T>
void SmartPointer<T>::allocate(int n)
{
    deallocate();

    n_ = n;

    // Check for 0 value, if 0 do not attempt an allocation

    if(n_ == 0)
        return;

    data_ = new T[n_];
}

template <class T>
void SmartPointer<T>::deallocate()
{
    if(data_ == NULL)
        return;

    delete[] data_;
    data_ = NULL;
    n_ = 0;
}

template <class T>
SmartPointer<T>& SmartPointer<T>::operator=(SmartPointer<T>& rhs)
{

    data_ = rhs.data_;
    n_ = rhs.n_;

    return *this;

}

template <class T>
SmartPointer<T>& SmartPointer<T>::operator=(T* rhs)
{

    data_ = rhs;
    n_ = 1;

    return *this;

}

template <class T>
SmartPointer<T>& SmartPointer<T>::operator=(T& rhs)
{

    data_ = &rhs;
    n_ = 1;

    return *this;

}

template <class T>
T* SmartPointer<T>::operator->()
{

    return data_;

}

template <class T>
T& SmartPointer<T>::operator()(int i)
{
    if(i < 0 || i >= n_)
        throw "Attempted to access element outside the bounds of SmartPointer";

    return data_[i];
}

#endif
