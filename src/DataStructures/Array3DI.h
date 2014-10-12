#ifndef ARRAY_3D_I_H
#define ARRAY_3D_I_H

#include <cstdlib>

#include "Array3D.h"

template <class T>
Array3D<T>::Array3D()
    :
      nI_(0),
      nJ_(0),
      nK_(0),
      n_(0),
      data_(NULL)
{

}

template <class T>
Array3D<T>::Array3D(int nI, int nJ, int nK)
    :
      Array3D()
{
    allocate(nI, nJ, nK);
}

template <class T>
Array3D<T>::~Array3D()
{
    deallocate();
}

template <class T>
void Array3D<T>::allocate(int nI, int nJ, int nK)
{
    int i, j;

    deallocate();

    nI_ = nI;
    nJ_ = nJ;
    nK_ = nK;
    n_ = nI_*nJ_*nK_;

    // Check for 0 value, if 0 do not attempt allocation

    if(n_ == 0)
        return;

    data_ = new T**[nI_];

    for(i = 0; i < nI_; ++i)
    {
        data_[i] = new T*[nJ_];
        
        for(j = 0; j < nJ_; ++j)
        {

            // Allocate nK_ + 1 so the iterator doesn't cause a seg fault
            if(i == 0 && j == 0)
            {

                data_[i][j] = new T[nK_];

            }
            else
            {

                data_[i][j] = new T[nK_ - 1];

            }
        }
    }
}

template <class T>
void Array3D<T>::deallocate()
{
    int i, j;

    if(data_ == NULL)
        return;

    for(i = 0; i < nI_; ++i)
    {
        for(j = 0; j < nJ_; ++j)
        {
            delete[] data_[i][j];
        }
        delete[] data_[i];
    }

    delete[] data_;

    data_ = NULL;

    nI_ = nJ_ = nK_ = n_ = 0;
}

template <class T>
T& Array3D<T>::operator()(int i, int j, int k)
{

    if(i < 0 || j < 0 || k < 0 ||
            i >= nI_ || j >= nJ_ || k >= nK_)
        throw "Attempted to access element outside the bounds of Array3D.";

    return data_[i][j][k];

}

//- Iterator methods

template <class T>
Array3D<T>::iterator::iterator()
    :
      dataPtr_(NULL),
      objectPtr_(NULL)
{

}

template <class T>
Array3D<T>::iterator::iterator(T* dataPtr,
                                      Array3D<T>* objectPtr,
                                      int i,
                                      int j,
                                      int k)
    :
      dataPtr_(dataPtr),
      objectPtr_(objectPtr),
      i_(i),
      j_(j),
      k_(k)
{

}

template <class T>
Array3D<T>::iterator& Array3D<T>::iterator::operator++()
{

    ++i_;

    if(i_ == objectPtr_->nI_)
    {
        i_ = 0;
        ++j_;
    }

    if(j_ == objectPtr_->nJ_)
    {
        j_ = 0;
        ++k_;
    }

    dataPtr_ = &objectPtr_->data_[i_][j_][k_];

    return *this;

}

template <class T>
T& Array3D<T>::iterator::operator*()
{

    return *dataPtr_;

}

template <class T>
bool Array3D<T>::iterator::operator!=(const iterator& rhs)
{

    if(dataPtr_ != rhs.dataPtr_)
        return true;

    return false;

}

template <class T>
Array3D<T>::iterator Array3D<T>::begin()
{

    return iterator(&data_[0][0][0], this, 0, 0, 0);

}

template <class T>
Array3D<T>::iterator Array3D<T>::end()
{

    // The end is at container size + 1 so that the iterator will iterate over the entire container

    return iterator(&data_[0][0][nK_], this, 0, 0, nK_);

}

#endif
