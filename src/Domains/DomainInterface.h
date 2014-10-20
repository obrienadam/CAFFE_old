#ifndef DOMAIN_INTERFACE_H
#define DOMAIN_INTERFACE_H

#include <vector>

#include "Input.h"
#include "SchemeInterface.h"

template <class STATE_TYPE>
class DomainInterface
{

protected:

public:

    virtual void initialize(Input& input) = 0;
    virtual int size() = 0;

    virtual void computeTimeDerivatives(STATE_TYPE* timeDerivatives) = 0;

};

#endif
