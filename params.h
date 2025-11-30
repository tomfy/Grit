// params.h  macros defining parameters controlling updates, etc.
// these are things we are likely to want to fiddle with
#define UPDATE_WHOLE_PATH FALSE
//#define NEW_MAKEREV TRUE

#define USE2EPSILONS FALSE
#define EPSILON_BIG_FACTOR 2.0 // Epsilon_big = Epsilon*EPSILON_BIG_FACTOR
#define BIGEPSILONFRACTION 0.5 // use Epsilon_big in update_path this fraction of updates. (If USE2EPSILONS == TRUE)

#define UPDATE_PATHS TRUE // chain.c run_chains.c

#define MAXSHORTPATHLENGTH 2 // up to this length, hits per path are kept track of.
 
// *****************************************
//#define USE_LOGRATIOS TRUE // TRUE -> work with logs for pi ratios, q ratios (to avoid overflow problems)

// lambda update macros
#define LAMBDA_DEF_FACTOR 1.0 // >= 1.0;
#define UPDATE_lambdaIT_MODE 0 // 0: lamdaI_prop = lambdaI + width*(drand() - 0.5), (reflected into 0 to lambda_max), else: use gamma_dev
#define LAMBDA_I_STEP_WIDTH_COEFF 6.0
#define LAMBDA_T_STEP_WIDTH_COEFF 18.0

#define rUPDATE_MODE 2 // 1: r, anything else: xi
#define r_STEP_WIDTH 100.0
#define XI_WIDTH 2.0

// seds preference macros
#define PREFERDD0SEDS TRUE //
//#define PREFERSEDSNDDM1MAX 0 // maximum number of delta_d=-1 choices for which prefer seds delta_d=0 (0 is recommended value jan10 03)
#define PDD0SEDS 0.5 // prob of doing a seds step when doing a delta_d = 0 step (if PREFERDD0SEDS and ..

#define PREFERDD0SEDS_SIGNFLIP TRUE
#define PDD0SEDS_SIGNFLIP 0.75 // prob of doing a seds step when doing a delta_d = 0 sign flip (if no delta_d = -1 steps possible).

// #define MAX_N_TEMPERATURES  10

#define PRINT_DD_INV_TRANS FALSE // Controls whether print how many delta_d = -1,0,1 inversions and translocations at each step along path

#define PREFER_INVERSIONS  TRUE // for delta d = 0, +1 steps
#define USE_IPF_LAMBDA_RATIO TRUE
#define IPF_EXPONENT 1.0
#define INITLONGPATH_IPF 10.0


#define COUNT_WITHIN_EQUIVALENCE_CLASS_STEPS FALSE // whether to histogram lengths which include steps which leave genome in same equivalence class or lengths which do not

#define LTHEAT FALSE  // TRUE to put in factor of LTHEATFACTOR^Lt for hottest chain in target distribution, for MCMCMC
#define LTHEATFACTOR 2.0 // see LTHEAT 
#define LTHEATKNEE 35 // up to this many translocations, enhance prob of high Lt paths in heated chains
#define LTSHARP  -1.0 // determines sharpness of cutoff of Lt heating near LTHEATKNEE
#define OLDLTHEAT FALSE
#define LTHEATBETA 1.5

#define FIXEDLAMBDAr FALSE
#define FIXED_LAMBDA 3.0 //  7123.36 //   578.0
#define FIXED_r 10000.0 // 400.0

#define NEWUPDATER TRUE
#define REV_OLDWAY FALSE
#define OLDWAYx FALSE

//#define OLDINIT TRUE
//#define RUNCHOLD FALSE // old or new version of run_chains

#define NSHORT (r_in->N_chain/2) // determines how many of the chains have a short initial path

#define SELF 1
#define ALL 2
#define NONE 3

#define CHROMENDMATCHING NONE // in counting "switches", chrom ends match: other chrom ends (SELF); everything (ALL), nothing (NONE)
//#define DDM1_PREFER_INVERSIONS TRUE  // FALSE to let each inv, trans be equally prob for delta_d = -1 steps
#define DDM1_IT_PROB_RATIO 2.0 // prob of choosing an inversion as next delta_d = -1 step (assuming there are both inv and trans to choose from)
//#define USE_DNSWITCHES TRUE
#define STEPFCN_USEDNSW FALSE
#define ACCEPTALL FALSE // forces acceptance of all proposed path updates, if TRUE

#define PROB_DDM1_PREFINV_USE_DNSWITCHES 0.3

#define SHORTCIRCUITCHROMSINCOMMON TRUE
#define DBDIST FALSE

#define UPDATE_TO_END FALSE
#define FLIPSIGNSTEPPROB 0.8 // this fraction of path update steps do a sign flip

//#define NEWGETABCD FALSE

// end  params.h
