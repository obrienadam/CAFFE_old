#ifndef FIELD_H
#define FIELD_H

#include <string>

#include "Array3D.h"

template <class T>
class Field : public Array3D<T>
{

public:

    Field(std::string name = "UnnamedField", bool isConserved = false)
        :
          name(name),
          isConserved(false)
    {

    }

    Field(int nI, int nJ, int nK, std::string name = "UnnamedField", bool isConserved = false)
        :
          Array3D<T>(nI, nJ, nK),
          name(name),
          isConserved(false)
    {

    }

    bool isConserved;
    std::string name;

};

#endif
