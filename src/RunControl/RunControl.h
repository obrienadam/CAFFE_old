#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include "ArgsList.h"

enum TerminationCondition{SIM_TIME, CPU_TIME, REAL_TIME, ITERATIONS};

class RunControl
{
  
 private:

  ArgsList argsList_;

  int itrs_;

  double simTime_, cpuTime_, realTime_;

 public:

  RunControl();
  RunControl(int argc, const char* argv[]);

  TerminationCondition terminationCondition;

  int maxItrs;
  
  double maxSimTime, maxCpuTime, maxRealTime;

  void displayStartMessage();

  bool continueRun(double timeStep = 0.);

};

#endif
