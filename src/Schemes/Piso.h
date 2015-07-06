#ifndef PISO_H
#define PISO_H

#include "Simple.h"

class Piso : public Simple
{
private:

    int nPCorrections_;

public:

    Piso();

    void initialize(Input &input, HexaFvmMesh &mesh);
    void discretize(double timeStep, std::vector<double>& timeDerivatives);
    void displayUpdateMessage();
};

#endif
