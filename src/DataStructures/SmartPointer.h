#ifndef SMART_POINTER_H
#define SMART_POINTER_H

#include <cstdlib>

template <class T>
class SmartPointer
{
  
 private:

  int n_;
  T* data_;

  bool isOriginalObject_;

 public:

  SmartPointer();
  SmartPointer(int n);
  SmartPointer(T* ptr);
  SmartPointer(const SmartPointer& other);
  ~SmartPointer();

  void allocate(int n);
  void deallocate();

  int size(){return n_;}

  void pushBack(const T& newData);

  SmartPointer& operator=(SmartPointer<T>& rhs);
  SmartPointer& operator=(T* rhs);
  SmartPointer& operator=(T& rhs);
  T* operator->();
  T& operator()(int);

};

#include "SmartPointerI.h"

#endif
