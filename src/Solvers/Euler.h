#ifndef EULER_H
#define EULER_H

#include "SolverInterface.h"

class Euler : public SolverInterface
{
  
 private:

 public:

  void advanceSolution(DomainInterface* domain, SchemeInterface* scheme);

};

#endif
