#ifndef TIME_H
#define TIME_H

#include <string>

#include <boost/date_time.hpp>

class Time
{
public:

    Time();
    Time(double hours, double minutes, double seconds);

    void tic();
    void toc();
    void pause();
    void resume();
    double elapsedTimeMicroseconds();
    std::string elapsedTime();
    std::string currentTime();
    bool exceededRefElapsedTime();

private:

    typedef boost::posix_time::ptime RealTime;
    typedef boost::posix_time::time_duration RealTimeDuration;

    RealTime ticTimePoint_, tocTimePoint_;
    RealTimeDuration elapsedTime_, refElapsedTime_;
};

#endif
