// for cpu time using vtimes (linux_gcc only)

// #if COMPILER == LINUX_GCC
// #include <sys/vtimes.h>
// #endif
#include <stdio.h>
#include <time.h>
#include <sys/types.h>

//#define NULL 0

double get_CPU_time(void);

double
get_CPU_time(void)
{
  fprintf(stderr, "CCCCCCCCCCCC: %g  %g\n", (double)clock(), (double)CLOCKS_PER_SEC);
  double x = (double)clock()/(double)CLOCKS_PER_SEC;
  fprintf(stderr, "time (get_CPU_time): %g\n", x);
  return x;
}
