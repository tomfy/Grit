// extdefs.h
// *************** external variable definitions ************************
// some external variables to keep track of memory allocation and freeing:

long n_step_alloc = 0;
long n_perm_alloc = 0;

long n_p_alloc = 0;
long n_path_alloc = 0;
long n_path_tree_node_alloc = 0;
long n_path_tree_alloc = 0;
long n_revseq_alloc = 0;

// other externals
Reversal r_init = {UNKNOWN,/*  UNKNOWN, */ UNKNOWN, UNKNOWN,
                   {UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN},
                   {UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN},
                   UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN, 
                   {UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN},
                   UNKNOWN, UNKNOWN, {UNKNOWN, UNKNOWN}}; 
                   
extern unsigned long rng_long; // defined in rngs.c, just declare extern here and in exts.h
int Error_count = 0; 
FILE* fp_in;
FILE* fp_out1, * fp_out2, * fp_out3, * fp_out4, * fp_out5, * fp_out6, * fp_out7, * fp_out8, * fp_out9, * fp_out10, * fp_rawhot, * fp_rawq;
FILE* fp_raw[MAX_N_TEMPERATURES]; 
FILE* fp_LiLtout, * fp_lambdaIlambdaTout;
int sigint_raised = 0; // flag; FALSE to start, TRUE after SIGINT signal handled
// ************* end of external variable definitions ***************
double run_around_clocks=0.0, run_around_calls=0.0;
double run_around_old_clocks=0.0, run_around_old_calls=0.0;
//double get_cycle_start_clocks=0.0,  get_cycle_start_calls=0.0;
double get_CD_clocks=0.0, get_CD_calls=0.0;
double get_CD_old_clocks=0.0, get_CD_old_calls=0.0;
int prop_lengths[2*MAXPATHLENGTH] = {0};
int acc[100][4] = {{0}};
int nstucks[100] = {0};
int lpa_hist[100] = {0};
int dd01_inv_all = 0, dd01_inv_acc = 0, dd01_inv_rej = 0;
int dd01_trans_all = 0, dd01_trans_acc = 0, dd01_trans_rej = 0;
int dd01_iort_all = 0, dd01_iort_acc = 0, dd01_iort_rej = 0;
int dd01_iandt_all = 0, dd01_iandt_acc = 0, dd01_iandt_rej = 0;

int n_path_p_a_nnoz = 0; // number of path updates with p_accept not normal or zero (i.e. is nan, inf, or denormal)
int n_lambdaIT_p_a_nnoz = 0; // number of lambdaI or lambdaT updates with log_p_accept not normal or zero (i.e. is nan, inf, or denormal)
int n_rxi_p_a_nnoz = 0;
int zeta_hist_a[20] = {0};
int zeta_hist_b[20] = {0};
long trans_a = 0, trans_b = 0;
int n_inv_prop = 0, n_trans_prop = 0;
int length_L_diff[50] = {0}; // length - L, where length includes steps that leave genome in same equivalence class, L does not include such
double ptl[MAX_N_MARKERS_PER_CHROMOSOME] = {0.0};
double atl[MAX_N_MARKERS_PER_CHROMOSOME] = {0.0};
//double dswfactor[5] = {1.0, 0.5, 0.25, 0.05, 0.025};
//double dswfactor[5] = {1.0, 1.0, 0.5, 0.05, 0.025};
double dswfactor[5] = {1.0, 0.7, 0.5, 0.2, 0.1};
//double dswfactor[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
int ddm1_prefer_inversions = TRUE, use_dnswitches = TRUE; // use to control kind of update step


//int dthhist[250] = {0};
int SMOOTHTYPE = 0;
//double smooth_width_over_stddev = 0.25;
