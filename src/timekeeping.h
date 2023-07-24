#ifndef _time_keeping_c_
#define _time_keeping_c_

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <time.h> //Already included in original LCM source code above
#include <sys/time.h>
#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */


// TIME MEASUREMENTS
#ifdef _WIN32
inline static double second()
{
	return ((double)clock()/(double)CLK_TCK);
}

#else

#include <sys/resource.h>

#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

inline static double myWallTime()
{
#ifdef __APPLE__
	static double timeConvert = 0.0;
	if ( timeConvert == 0.0 )
	{
		mach_timebase_info_data_t timeBase;
		mach_timebase_info(&timeBase);
		timeConvert = (double)timeBase.numer / (double)timeBase.denom / 1000000000.0;
	}
	return mach_absolute_time() * timeConvert;
#else
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + 1.0e-9*((double)ts.tv_nsec);
#endif // __APPLE__
}

inline static double second()
{
	double t = myWallTime();
	return(t);
}

#endif
// END TIME MEASUREMENTS



// Measure running time
double measureTime(){
  struct rusage t;
  struct timeval tv,ts;
  getrusage(RUSAGE_SELF, &t);
  tv = t.ru_utime;
  ts = t.ru_stime;
  return tv.tv_sec + ts.tv_sec + ((double)tv.tv_usec + (double)ts.tv_usec) * 1e-6;
}

// Measure peak memory usage
size_t measurePeakMemory(){
  struct rusage t;
  getrusage(RUSAGE_SELF, &t);
  return (size_t)t.ru_maxrss;
}



#endif
