#include "timeBasedSeed.hh"
#include <sys/time.h>

unsigned int timeBasedSeed() {

  static int previous_value = 0; 
  static const int repeat_preventer = 171;
  timeval a;
  gettimeofday(&a,0); 
  int sec = a.tv_sec;
  int usec = a.tv_usec;
  sec = (sec%4000) - 2000;
  usec = usec%1000000;
  int retval = 1000000*sec+usec;
  if ( retval == previous_value) retval += repeat_preventer;
  previous_value = retval;
  return sec+usec;

}
