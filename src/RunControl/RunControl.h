#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include <string>

#include <boost/date_time.hpp>
#include <boost/chrono.hpp>

#include "ArgsList.h"
#include "Input.h"

#include "DomainIncludes.h"
#include "SolverIncludes.h"
#include "SchemeIncludes.h"

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

    //- Set the run control parameters

    void setRunControlParametersFromInputFile();

    //- A function that allows an external check for vaild command line arguments

    bool validOptionsSelected();

    //- Evaluate whether or not run should continue

    bool continueRun(double timeStep = 0.);

    //- Output messages

    void displayStartMessage();
    void displayUpdateMessage();
    void displayEndMessage();

    //- Initializes the solver to a specific type

    template <class DOMAIN_TYPE, class STATE_TYPE>
    void solverInitialize(SolverInterface<DOMAIN_TYPE, STATE_TYPE>*& solver);

};

#include "RunControlI.h"

#endif
