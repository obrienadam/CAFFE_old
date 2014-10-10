#ifndef EULER_H
#define EULER_H

#include "SolverInterface.h"

class Euler : public SolverInterface
{
  
 private:

 public:

  Euler();

  void advanceSolution(DomainInterface* domain, SchemeInterface* scheme);

  void initialize(Input &input);

};

#endif
