#ifndef SMART_POINTER_H
#define SMART_POINTER_H

template <class T>
class SmartPointer
{
  
 private:

  int n_;
  T* data_;

 public:

  SmartPointer();
  SmartPointer(int n);
  ~SmartPointer();

  void allocate(int n);
  void deallocate();

  int size(){return n_;}

  T& operator()(int i);

};

#include "SmartPointerI.h"

#endif
