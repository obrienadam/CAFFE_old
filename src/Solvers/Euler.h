#ifndef EULER_H
#define EULER_H

#include "SolverInterface.h"

class Euler : public SolverInterface
{

private:

public:

    Euler();

    void initialize(Input &input);

    void advanceSolution(SmartPointer<DomainInterface> domain, SmartPointer<SchemeInterface> scheme);

};

#endif
