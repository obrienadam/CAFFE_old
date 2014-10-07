#include <cstdlib>

template <class T>
SmartPointer3D<T>::SmartPointer3D() 
:
nI_(0),
  nJ_(0),
  nK_(0),
  n_(0),
  data_(NULL)
{

}

template <class T>
SmartPointer3D<T>::SmartPointer3D(int nI, int nJ, int nK)
:
SmartPointer3D()
{
  allocate(nI, nJ, nK);
}

template <class T>
SmartPointer3D<T>::~SmartPointer3D()
{
  deallocate();
}

template <class T>
void SmartPointer3D<T>::allocate(int nI, int nJ, int nK)
{
  int i, j;

  deallocate();

  nI_ = nI;
  nJ_ = nJ;
  nK_ = nK;
  n_ = nI_*nJ_*nK_;

  data_ = new T**[nI_];

  for(i = 0; i < nI_; ++i)
    {
      data_[i] = new T*[nJ_];
        
      for(j = 0; j < nJ_; ++j)
        {
          data_[i][j] = new T[nK_];
        }
    }
}

template <class T>
void SmartPointer3D<T>::deallocate()
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
T& SmartPointer3D<T>::operator()(int i, int j, int k)
{
  if(i < 0 || j < 0 || k < 0 ||
     i >= nI_ || j >= nJ_ || k >= nK_)
    throw "Attempted to access element outside the bounds of SmartPointer3D.";

  return data_[i][j][k];
}
