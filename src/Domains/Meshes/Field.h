#ifndef FIELD_H
#define FIELD_H

#include <string>

#include "Array3D.h"

enum{CONSERVED, AUXILLARY};

template <class T>
class Field : public Array3D<T>
{

public:

    Field(std::string name = "UnnamedField", int type = AUXILLARY)
        :
          name(name),
          type(type)
    {

    }

    Field(int nI, int nJ, int nK, std::string name = "UnnamedField", int type = AUXILLARY)
        :
          Array3D<T>(nI, nJ, nK),
          name(name),
          type(type)
    {

    }

    //- The "type" determines whether or not tranport equations need to be solved for this field

    bool type;
    std::string name;

};

#endif
