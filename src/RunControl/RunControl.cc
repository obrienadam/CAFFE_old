#include <iostream>

#include "RunControl.h"

RunControl::RunControl()
 :
  itrs_(0),
  simTime_(0.),
  cpuTime_(0.),
  realTime_(0.)
{

}

RunControl::RunControl(int argc, const char* argv[])
:
  RunControl()
{
  argsList_.readArgs(argc, argv);
}

bool RunControl::continueRun(double timeStep)
{
  ++itrs_;
  simTime_ += timeStep;

  switch(terminationCondition)
    {

    case ITERATIONS:

      if(itrs_ >= maxItrs)
	return false;

      break;


    case SIM_TIME:
      
      if(simTime_ >= maxSimTime)
	return false;

      break;

    case CPU_TIME:

      if(cpuTime_ >= maxCpuTime)
	return false;

      break;

    case REAL_TIME:

      if(realTime_ >= maxRealTime)
	return false;

      break;

    default:

      throw "Invalid termination condition in RunControl.";
      return false;

    };

      return true;
}

void RunControl::displayStartMessage()
{
  using namespace std;

  cout << endl 
       << "#################################################################\n"
       << endl
       << "|  || ___\n" 
       << "|  ||/ _  \\ \\\n"				\
       << "|  ||\\__/ | | |\n"
       << "|  | \\___/  | | CAFFE\n"
       << "|   \\______/  /\n"
       << " \\___________/\n"  
       << endl
       << "  Computational Algorithm Framework for Fluid Equations (CAFFE)\n"  
       << endl
       << "                      Author: Adam O'Brien\n"
       << "	            E-mail: roni511@gmail.com\n"
       << endl
       << "#################################################################\n";
  
}
