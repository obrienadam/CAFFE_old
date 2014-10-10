#ifndef SCHEME_INTERFACE_H
#define SCHEME_INTERFACE_H

#include "Input.h"

class SchemeInterface
{

 private:

 public:

    virtual void initialize(Input& input) = 0;

};

#endif
