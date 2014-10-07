#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include <boost/date_time.hpp>
#include <boost/chrono.hpp>

#include "ArgsList.h"
#include "Input.h"

enum TerminationCondition{SIM_TIME, CPU_TIME, REAL_TIME, ITERATIONS};

class RunControl
{

 typedef boost::posix_time::ptime RealTime;
 typedef boost::posix_time::time_duration RealTimeDuration;
 typedef boost::chrono::process_cpu_clock::duration CpuTimeDuration;

 private:

  ArgsList argsList_;
  Input input_;

  int itrs_;

  double simTime_, maxSimTime_;

  RealTime startTime_;
  RealTimeDuration elapsedTime_, maxElapsedTime_;
  CpuTimeDuration cpuTime_, maxCpuTime_;

 public:

  RunControl();
  RunControl(int argc, const char* argv[]);

  TerminationCondition terminationCondition;

  int maxItrs;

  void displayStartMessage();

  bool continueRun(double timeStep = 0.);

  void displayUpdateMessage();

  void displayEndMessage();

};

#endif
