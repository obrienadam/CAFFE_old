#include <sstream>

#include "RunControl.h"
#include "Output.h"

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
        Output::raiseException("RunControl", "continueRun", "invalid termination condition \"" + terminationCondition_ + "\" selected.");
    }

    return true;
}

void RunControl::displayStartMessage()
{
    using namespace std;
    using namespace boost::posix_time;

    ostringstream message;

    startTime_ = second_clock::local_time();

    Output::printLine();

    message << "Beginning simulation. Terminating on condition: " << terminationCondition_ << ".";

    Output::printToScreen(message.str());

    message.str("");

    message << "Iterations beginning on " << startTime_ << ".";

    Output::printToScreen(message.str());
}

void RunControl::displayUpdateMessage()
{
    std::ostringstream message;
    double completionPercentage;

    if(terminationCondition_ == "iterations")
    {
        completionPercentage = double(itrs_)/double(maxItrs_)*100.;
    }
    else if(terminationCondition_ == "simTime")
    {
        completionPercentage = simTime_/maxSimTime_*100.;
    }
    else if(terminationCondition_ == "cpuTime")
    {

    }
    else if (terminationCondition_ == "realTime")
    {

    }

    message << "Simulation completion: " << completionPercentage << "%\n"
            << "Iterations completed: " << itrs_ << "\n"
            << "Simulation time: " << simTime_ << "\n"
            << "Elapsed time: " << elapsedTime_ << "\n"
            << "CPU time: " << cpuTime_;

    Output::printToScreen(message.str());
}

void RunControl::displayEndMessage()
{
    using namespace std;

    ostringstream message;

    Output::printLine();

    message << "Iterations complete on " << startTime_ + elapsedTime_;

    Output::printToScreen(message.str());

    message.str("");

    message << "Iterations completed: " << itrs_;

    Output::printToScreen(message.str());

    message.str("");

    message << "Simulation time: " << simTime_;

    Output::printToScreen(message.str());

    message.str("");

    message << "Elapsed time: " << elapsedTime_;

    Output::printToScreen(message.str());

    message.str("");

    message << "CPU time: " << cpuTime_;

    Output::printToScreen(message.str());
    Output::printLine();
}

void RunControl::initializeCase(Solver& solver, DomainInterface& domain)
{
    domain.initialize(input_);
    solver.initialize(input_);
}
