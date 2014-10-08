#include <iostream>

#include "RunControl.h"

RunControl::RunControl()
    :
      itrs_(0),
      simTime_(0.),
      startTime_(boost::posix_time::microsec_clock::local_time()),
      elapsedTime_(0, 0, 0, 0),
      maxElapsedTime_(48, 0, 0, 0),
      terminationCondition_("iterations")
{

}

RunControl::RunControl(int argc, const char* argv[])
    :
      RunControl()
{

    argsList_.readArgs(argc, argv);
    input_.openInputFile(argsList_.inputFilename_);
    setRunControlParametersFromInputFile();

}

void RunControl::setRunControlParametersFromInputFile()
{

    terminationCondition_ = input_.inputStrings["terminationCondition"];
    maxItrs_ = input_.inputInts["maxItrs"];
    maxSimTime_ = input_.inputDoubles["maxSimTime"];

}

bool RunControl::continueRun(double timeStep)
{

    ++itrs_;
    simTime_ += timeStep;

    if(terminationCondition_ == "iterations")
    {

        if(itrs_ >= maxItrs_)
            return false;

    }

    else if(terminationCondition_ == "simTime")
    {

        if(simTime_ >= maxSimTime_)
            return false;

    }

    else if(terminationCondition_ == "cpuTime")
    {

        if(cpuTime_ >= maxCpuTime_)
            return false;

    }

    else if (terminationCondition_ == "realTime")
    {

        if(elapsedTime_ >= maxElapsedTime_)
            return false;

    }

    else
    {

        std::string errorMessage("Invalid termination condition \"" + terminationCondition_ + "\" selected in RunControl.");

        throw errorMessage.c_str();

    }

    return true;

}

void RunControl::displayStartMessage()
{

    using namespace std;
    using namespace boost::posix_time;

    startTime_ = second_clock::local_time();

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
         << "Iterations beginning on " << startTime_ << ".\n"
         << endl
         << "Simulation termination condition: " << terminationCondition_ << endl
         << endl;

}

void RunControl::displayUpdateMessage()
{

}

void RunControl::displayEndMessage()
{

    using namespace std;

    cout << "Iterations complete on " << startTime_ + elapsedTime_ << ".\n"
         << endl
         << "Elapsed time: " << elapsedTime_ << endl
         << "CPU time: " << cpuTime_ << endl;

}
