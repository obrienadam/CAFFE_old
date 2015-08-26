#include <cstdlib>

#include "Array3D.h"

template <class T>
Array3D<T>::Array3D(int sizeI, int sizeJ, int sizeK)
    :
      std::vector<T>::vector(sizeI*sizeJ*sizeK)
{
    setSizes(sizeI, sizeJ, sizeK);
}

template <class T>
void Array3D<T>::resize(int sizeI, int sizeJ, int sizeK)
{
    setSizes(sizeI, sizeJ, sizeK);
    std::vector<T>::resize(sizeIJ_*sizeK_);
}

template <class T>
void Array3D<T>::clear()
{
    std::vector<T>::clear();
    sizeI_ = sizeJ_ = sizeK_ = sizeIJ_ = 0;
}

template<class T>
void Array3D<T>::assign(const T &val)
{
    std::vector<T>::assign(std::vector<T>::size(), val);
}

template<class T>
void Array3D<T>::add(const T &val)
{
    for(auto it = std::vector<T>::begin(); it != std::vector<T>::end(); ++it)
        *it += val;
}

template <class T>
T& Array3D<T>::operator()(int i, int j, int k)
{
    if(i < 0 || j < 0 || k < 0 ||
            i >= sizeI_ || j >= sizeJ_ || k >= sizeK_)
        throw "Attempted to access element outside the bounds of Array3D.";

    return std::vector<T>::operator [](k*sizeIJ_ + j*sizeI_ + i);
}

template <class T>
const T& Array3D<T>::operator()(int i, int j, int k) const
{
    if(i < 0 || j < 0 || k < 0 ||
            i >= sizeI_ || j >= sizeJ_ || k >= sizeK_)
        throw "Attempted to access element outside the bounds of Array3D.";

    return std::vector<T>::operator [](k*sizeIJ_ + j*sizeI_ + i);
}

/***************************** Protected methods ********************************/

template <class T>
void Array3D<T>::setSizes(int sizeI, int sizeJ, int sizeK)
{
    sizeI_ = sizeI;
    sizeJ_ = sizeJ;
    sizeK_ = sizeK;
    sizeIJ_ = sizeI_*sizeJ_;
}
