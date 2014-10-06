#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <string>

#include "DomainInterface.h"
#include "SchemeInterface.h"

class SolverInterface
{

 protected:

  std::string solverName_;

 public:

  virtual void advanceSolution(DomainInterface* domain, SchemeInterface* scheme) = 0;

};

#endif
