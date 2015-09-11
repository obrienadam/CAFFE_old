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

    double elapsedTimeMicroseconds() const;
    std::string elapsedTime() const;
    std::string currentTime() const;
    bool exceededRefElapsedTime() const;

private:

    typedef boost::posix_time::ptime RealTime;
    typedef boost::posix_time::time_duration RealTimeDuration;

    RealTime ticTimePoint_, tocTimePoint_;
    RealTimeDuration elapsedTime_, refElapsedTime_;
};

#endif
