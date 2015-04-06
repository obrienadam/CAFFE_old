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
    deallocate();

    nI_ = nI;
    nJ_ = nJ;
    nK_ = nK;
    nInJ_ = nI_*nJ_;
    nInK_ = nI_*nK_;
    nJnK_ = nJ_*nK_;
    n_ = nInJ_*nK_;

    // Check for 0 value, if 0 do not attempt allocation

    if(n_ == 0)
        return;
    
    data_ = new T[n_ + 1];
}

template <class T>
void Array3D<T>::deallocate()
{
    if(data_ == NULL)
        return;

    delete[] data_;
    data_ = NULL;

    nI_ = nJ_ = nK_ = nInJ_ = nInK_ = nJnK_ = n_ = 0;
}

template <class T>
T& Array3D<T>::operator()(int i, int j, int k)
{
    if(i < 0 || j < 0 || k < 0 ||
            i >= nI_ || j >= nJ_ || k >= nK_)
        throw "Attempted to access element outside the bounds of Array3D.";

    return data_[k*nInJ_ + j*nI_ + i];
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
                                      int k)
    :
      dataPtr_(dataPtr),
      objectPtr_(objectPtr),
      k_(k)
{

}

template <class T>
typename
Array3D<T>::iterator& Array3D<T>::iterator::operator++()
{
    ++k_;

    dataPtr_ = &objectPtr_->data_[k_];

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
typename
Array3D<T>::iterator Array3D<T>::begin()
{
    return iterator(&data_[0], this, 0);
}

template <class T>
typename
Array3D<T>::iterator Array3D<T>::end()
{
    // The end is at container size + 1 so that the iterator will iterate over the entire container

    return iterator(&data_[n_], this, n_);
}
