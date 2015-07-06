#ifndef PISO_H
#define PISO_H

#include "Simple.h"

class Piso : public Simple
{
public:

    void discretize(double timeStep, std::vector<double>& timeDerivatives);

    void displayUpdateMessage();
};

#endif
