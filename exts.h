// exts.h

// declarations of external variables
// some external variables to keep track of memory allocation and freeing:

extern long n_step_alloc;
extern long n_perm_alloc;
extern long n_p_alloc;
extern long n_path_alloc;
extern long n_path_tree_node_alloc ;
extern long n_path_tree_alloc ;
extern long n_revseq_alloc ;

// other externals
extern Reversal r_init;
extern unsigned long rng_long;
extern int Error_count;
extern FILE* fp_in;
extern FILE* fp_out1, * fp_out2, * fp_out3, * fp_out4, * fp_out5,
  * fp_out6, * fp_out7, * fp_out8, * fp_rawhot, * fp_rawq, * fp_prog;
extern FILE* fp_raw[MAX_N_TEMPERATURES]; 
extern FILE* fp_LiLtout, * fp_lambdaIlambdaTout;
extern int sigint_raised; // flag; FALSE to start, TRUE after SIGINT signal handled
extern char* output_prefix;

extern double run_around_clocks, run_around_calls;
extern double run_around_old_clocks, run_around_old_calls;
extern double get_CD_clocks, get_CD_calls;
extern double get_CD_old_clocks, get_CD_old_calls;
extern int prop_lengths[2*MAXPATHLENGTH];
extern int acc[100][4];
extern int nstucks[100];
extern int lpa_hist[100];
extern int dd01_inv_all, dd01_inv_acc, dd01_inv_rej;
extern int dd01_trans_all, dd01_trans_acc, dd01_trans_rej;
extern int dd01_iort_all, dd01_iort_acc, dd01_iort_rej;
extern int dd01_iandt_all, dd01_iandt_acc, dd01_iandt_rej;

extern int n_path_p_a_nnoz; // number of path updates with p_accept not normal or zero (i.e. is nan, inf, or denormal)
extern int n_lambdaIT_p_a_nnoz; // number of lambdaI or lambdaT updates with p_accept not normal or zero (i.e. is nan, inf, or denormal)
extern int n_rxi_p_a_nnoz;
extern int zeta_hist_a[20];
extern int zeta_hist_b[20];
extern long trans_a, trans_b;
extern int n_inv_prop, n_trans_prop;
extern int length_L_diff[40]; //
extern double ptl[MAX_N_MARKERS_PER_CHROMOSOME];
extern double atl[MAX_N_MARKERS_PER_CHROMOSOME];
extern double dswfactor[5];
extern int ddm1_prefer_inversions, use_dnswitches; // use to control kind of update step

extern int SMOOTHTYPE;
//extern double smooth_width_over_stddev;

//extern int dthhist[250];
/* extern long delta_m_req[3]; */
/* extern long delta_ci_prop[3]; */
/* extern long dmci[3][3]; */

/* extern long dci_nb_m1[3][7]; */
/* extern long dci_nb_0[3][7]; */
/* extern long dci_nb_1[3][7]; */

// extern int use_cicg; // if == TRUE, use c_i and c_g in generating paths, etc. else use c_i and m_ne

// end of exts.h
