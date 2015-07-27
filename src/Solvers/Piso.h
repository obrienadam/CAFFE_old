#ifndef PISO_H
#define PISO_H

#include "Simple.h"

class Piso : public Simple
{
public:

    Piso(const Input &input, const HexaFvmMesh &mesh);

    virtual double solve(double timeStep);
    virtual void displayUpdateMessage();

protected:

    int nPCorrections_;
};

#endif
