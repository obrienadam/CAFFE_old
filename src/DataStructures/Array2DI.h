#include <cstdlib>

template <class T>
Array2D<T>::Array2D()
    :
      nI_(0),
      nJ_(0),
      n_(0),
      data_(NULL)
{

}

template <class T>
Array2D<T>::Array2D(int nI, int nJ)
    :
      Array2D()
{
    allocate(nI, nJ);
}

template <class T>
Array2D<T>::~Array2D()
{
    deallocate();
}

template <class T>
void Array2D<T>::allocate(int nI, int nJ)
{
    deallocate();

    nI_ = nI;
    nJ_ = nJ;
    n_ = nI_*nJ_;

    // Check for 0 value, if 0 do not attempt allocation

    if(n_ == 0)
        return;
    
    data_ = new T[n_ + 1];
}

template <class T>
void Array2D<T>::deallocate()
{
    if(data_ == NULL)
        return;

    delete[] data_;
    data_ = NULL;

    nI_ = nJ_ = n_ = 0;
}

template <class T>
T& Array2D<T>::operator()(int i, int j)
{
    if(i < 0 || j < 0 ||
            i >= nI_ || j >= nJ_)
        throw "Attempted to access element outside the bounds of Array2D.";

    return data_[j*nI_ + i];
}

//- Iterator methods

template <class T>
Array2D<T>::iterator::iterator()
    :
      dataPtr_(NULL),
      objectPtr_(NULL)
{

}
