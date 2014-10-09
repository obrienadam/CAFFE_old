#ifndef SMART_POINTER_H
#define SMART_POINTER_H

#include <cstdlib>

#include "SmartPointer.h"

template <class T>
SmartPointer<T>::SmartPointer()
    :
      n_(0),
      data_(NULL)
{

}

template <class T>
SmartPointer<T>::SmartPointer(int n)
    :
      data_(NULL)
{
    allocate(n);
}

template <class T>
SmartPointer<T>::~SmartPointer()
{
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
T& SmartPointer<T>::operator()(int i)
{
    if(i < 0 || i >= n_)
        throw "Attempted to access element outside the bounds of SmartPointer";

    return data_[i];
}

#endif
