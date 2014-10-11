#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <string>
#include <vector>

#include "DomainInterface.h"
#include "SchemeInterface.h"
#include "Input.h"

class SolverInterface
{

protected:

    SolverInterface(std::string solverName)
        :
          solverName_(solverName)
    {

    }

    std::string solverName_;

    double timeStep_;

    std::vector<DomainInterface> timeDerivatives_;

public:

    virtual void advanceSolution(DomainInterface* domain, SchemeInterface* scheme) = 0;

    virtual void initialize(Input& input) = 0;

};

#endif
