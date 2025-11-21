// grit.h


#define NDEBUG  // comment out to turn on assert
#define MYDEBUG

#define DOTESTHITS FALSE
#define TRANSITION_TESTS FALSE
#define ERROR_COUNT_LIMIT 1
#define TRUE 1
#define FALSE 0
#define RIGHT 0
#define LEFT 1
#define UNKNOWN -10000 //set cycle_number, etc. to this to indicate not calculated yet
#define RCFACTOR_1 200.0
#define RCFACTOR_2 5000.0
//#define RNG 5 // choose random number generator; now done in rng.c

#define LEFTTORIGHT TRUE
//#define N_SIGNED_TO_TRY 100 // to get unsigned inversion distance, do this many signed perms (get_signs_short)

// #define NOTSTOPPROB (4.0*Epsilon*Epsilon)  // moved to path.c
// #define BETA 0.0 // moved to path.c
// #define EPSILON2FACTOR 0.5 // moved to perm.c

#define WRITERAW TRUE
#define WRITERAWEVERY 8 // write raw data this often (every this many path updates)
#define MYFPNORMAL 4  // Determined empirically. I don't know why FP_NORMAL is not defined (defined in math.h)
#define HANDLE_SIGINT TRUE

#define EITHER 1     // EITHER, INVERSION, TRANSLOCATION 
#define INVERSION 2
#define TRANSLOCATION 3

#define MAX_N_CHROMOSOMES 40
#define MAX_N_MARKERS_PER_CHROMOSOME 500
#define MAXMARKERNAMELENGTH 16

// run parameter limits
#define NMARKERMIN 4
#define NMARKERMAX 10000
#define NCHAINMIN 1
#define NCHAINMAX 80
#define NDATAMIN 1
#define NDATAMAX 1000 
#define MINRUNLENGTH 10
#define MAXRUNLENGTH 100000000
#define MAX_N_TEMPERATURES 10

#define MAXPATHLENGTH 600

// output related macros
// avg (L etc.) over Navg steps and output every Navg MC steps
#define NAVG0 16 // initially Navg = NAVG0
#define N_DOUBLE 64  // every N_DOUBLE outputs, Navg gets doubled


#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <malloc.h> // 

#define LINUX_GCC 1
#define MAC_GCC 2
#define WINDOWS_BORLAND 3
#define COMPILER LINUX_GCC

#if COMPILER == LINUX_GCC

#include <stdio.h>
#include <time.h>
#include <float.h>
//#include <sys/vtimes.h>   
#define RANDMAX RAND_MAX
#define RANDF rand
#define SEEDRAND srand
#define WAITFORKEY ;   //  do nothing
#define GETCPUTIME  (get_CPU_time())
// double get_CPU_time();
#define SM_FACT 10.0 // number to multiply DBL_EPSILON ( = 2.22e-16) by when checking that something is small)

#elif COMPILER == MAC_GCC // MAC OSX GCC (BEAZY'S COMPUTER)

#include <stdio.h>
#include <time.h>
#include <float.h>
//#include <sys/vtimes.h>   link get_CPU_time.o
#define RANDMAX RAND_MAX
#define RANDF rand
#define WAITFORKEY ;   //  do nothing
#define SEEDRAND srand#
#define GETCPUTIME (0.0)
//double get_CPU_time();
#define SM_FACT 10.0 // number to multiply DBL_EPSILON ( = 2.22e-16) by when checking that something is small)

#else // borland

#include <stdio.h>
#include <time.h>
    //#define RANDMAX LRAND_MAX
    //#define RANDF _lrand
    //#define SEEDRAND srand
#define WAITFORKEY getchar();  // wait for enter key (to keep dos window from disappearing
#define GETCPUTIME 0.0
#define DBL_EPSILON 2.2e-16

#endif

#define SHORTINT int

// because MS c++ doesn't define M_PI (why?)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void check_alloc_info(FILE* fp);
void handle_control_c(int sig);
// void process_error(void); // counts errors, exits if too many
int compare_ints(const void* aa, const void* bb);
int compare_doubles(const void* aa, const void* bb);
double trimean(const int* array, int n);
double rmean(int* Nconv, int n, double* sigma, double f);
void print_other_output(void);

// end of grit.h
