// for cpu time using vtimes (linux_gcc only)

// #if COMPILER == LINUX_GCC
// #include <sys/vtimes.h>
// #endif

#define NULL 0

double get_CPU_time();

double
get_CPU_time()
{
     return 0.0;
    
/*    struct vtimes vt;
    vtimes(&vt, NULL);
    return (double)(vt.vm_utime + vt.vm_stime)/(double)VTIMES_UNITS_PER_SECOND;
*/
    }
