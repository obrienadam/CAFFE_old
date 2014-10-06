#ifndef SMART_POINTER_3D_H
#define SMART_POINTER_3D_H

template <class T>
class SmartPointer3D
{

 private:

  int nI_, nJ_, nK_, n_;
  T*** data_;

 public:

  SmartPointer3D();
  SmartPointer3D(int nI, int nJ, int nK);
  ~SmartPointer3D();

  void allocate(int nI, int nJ, int nK);
  void deallocate();

  int sizeI(){return nI_;}
  int sizeJ(){return nJ_;}
  int sizeK(){return nK_;}
  int size(){return n_;}
  
  T& operator()(int i, int j, int k);
};

#include "SmartPointer3DI.h"

#endif
