#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include <string>

#include <boost/date_time.hpp>
#include <boost/chrono.hpp>

#include "ArgsList.h"
#include "Input.h"

#include "DomainInterface.h"
#include "Solver.h"

class RunControl
{
    //- Typedefs

    typedef boost::posix_time::ptime RealTime;
    typedef boost::posix_time::time_duration RealTimeDuration;
    typedef boost::chrono::process_cpu_clock::duration CpuTimeDuration;

private:

    //- Input

    ArgsList argsList_;
    Input input_;

    //- Simulation control

    std::string terminationCondition_;
    int itrs_, maxItrs_;
    double simTime_, maxSimTime_;

    //- Time related objects

    RealTime startTime_;
    RealTimeDuration elapsedTime_, maxElapsedTime_;
    CpuTimeDuration cpuTime_, maxCpuTime_;

    //- Private constructor only used to initialize defaults

    RunControl();

public:

    //- Constructor that accepts command line arguments

    RunControl(int argc, const char* argv[]);

    //- Evaluate whether or not run should continue

    bool continueRun(double timeStep = 0.);

    //- Output messages

    void displayStartMessage();
    void displayUpdateMessage();
    void displayEndMessage();

    //- Initializes the various case objects

    void initializeCase(Solver& solver, DomainInterface& domain);
};

#endif
