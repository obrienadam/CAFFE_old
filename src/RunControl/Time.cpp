#include "Time.h"

Time::Time()
{
    tic();
}

Time::Time(double hours, double minutes, double seconds)
    :
      refElapsedTime_(hours, minutes, seconds, 0.)
{

}

void Time::tic()
{
    ticTimePoint_ = boost::posix_time::microsec_clock::local_time();
    elapsedTime_ = RealTimeDuration(0., 0., 0., 0.);
}

void Time::toc()
{
    tocTimePoint_ = boost::posix_time::microsec_clock::local_time();
    elapsedTime_ += tocTimePoint_ - ticTimePoint_;
}

void Time::pause()
{
    elapsedTime_ += tocTimePoint_ - ticTimePoint_;
    ticTimePoint_ = boost::posix_time::microsec_clock::local_time();
    tocTimePoint_ = boost::posix_time::microsec_clock::local_time();
}

void Time::resume()
{
    ticTimePoint_ = boost::posix_time::microsec_clock::local_time();
}

double Time::elapsedTimeMicroseconds() const
{
    return elapsedTime_.total_microseconds();
}

std::string Time::elapsedTime() const
{
    return boost::posix_time::to_simple_string(elapsedTime_);
}

std::string Time::currentTime() const
{
    return boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time());
}

bool Time::exceededRefElapsedTime() const
{
    return elapsedTime_ > refElapsedTime_;
}

