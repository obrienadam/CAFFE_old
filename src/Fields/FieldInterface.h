#ifndef FIELD_INTERFACE_H
#define FIELD_INTERFACE_H

#include "SmartPointer3D.h"
#include "Point3D.h"

template <class T>
class FieldInterface
{

 private:

  SmartPointer3D<T> field_;
  SmartPointer3D<Point3D> coords_;

 public:

  FieldInterface(){}
 FieldInterface(int nI, int nJ, int nK)
   :
  field_(nI, nJ, nK),
    coords_(nI, nJ, nK)
      {

      }

  T& operator()(int i, int j, int k)
    {
      return field_(i, j, k);
    }

};

#endif
