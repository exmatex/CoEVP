#include <ctime>

#include "fastTimer.hpp"

double getUnixTime(void)
{
    struct timespec ts;

    if(clock_gettime(CLOCK_REALTIME, &ts) != 0) return 0;

    return (((double) ts.tv_sec) + (double) (ts.tv_nsec / 1000000000.0));
}

