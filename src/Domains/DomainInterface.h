#ifndef DOMAIN_INTERFACE_H
#define DOMAIN_INTERFACE_H

#include <vector>

#include "Input.h"
#include "SchemeInterface.h"


class DomainInterface
{

protected:

public:

    virtual void allocate(Input& input) = 0;
    virtual int size() = 0;



};

#endif
