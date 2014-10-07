#include <iostream>

#include "RunControl.h"

RunControl::RunControl()
 :
  itrs_(0),
  simTime_(0.),
  startTime_(boost::posix_time::microsec_clock::local_time()),
  elapsedTime_(0, 0, 0, 0),
  maxElapsedTime_(48, 0, 0, 0),
  terminationCondition(ITERATIONS)
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
      
      if(simTime_ >= maxSimTime_)
	return false;

      break;

    case CPU_TIME:

      if(cpuTime_ >= maxCpuTime_)
	return false;

      break;

    case REAL_TIME:

      if(elapsedTime_ >= maxElapsedTime_)
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

  startTime_ = boost::posix_time::second_clock::local_time();

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
       << "#################################################################\n"
       << endl
       << "Iterations beginning on " << startTime_ << "." << endl
       << endl;
  
}

void RunControl::displayUpdateMessage()
{

}

void RunControl::displayCompletionMessage()
{

}
