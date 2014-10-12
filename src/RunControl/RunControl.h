#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include <string>

#include <boost/date_time.hpp>
#include <boost/chrono.hpp>

#include "ArgsList.h"
#include "Input.h"
#include "SmartPointer.h"

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

    ArgsList argsList_;
    Input input_;

    std::string terminationCondition_;

    int itrs_, maxItrs_;

    double simTime_, maxSimTime_;

    RealTime startTime_;
    RealTimeDuration elapsedTime_, maxElapsedTime_;
    CpuTimeDuration cpuTime_, maxCpuTime_;

    RunControl();

public:

    RunControl(int argc, const char* argv[]);

    void setRunControlParametersFromInputFile();

    bool validOptionsSelected();

    bool continueRun(double timeStep = 0.);

    //- Output messages

    void displayStartMessage();
    void displayUpdateMessage();
    void displayEndMessage();

    //- Initialization

    void initializeObjects(SmartPointer<DomainInterface> domain,
                           SmartPointer<SolverInterface> solver,
                           SmartPointer<SchemeInterface> scheme);

};

#endif
