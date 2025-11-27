// structs.h    structure definitions for grit

typedef struct real_pair{
    double a, b;
}Real_pair;

typedef struct int_pair{
    int a, b;
}Int_pair;

typedef struct four_longs{
    long a, b, c, d;
}Four_longs;

typedef struct four_doubles{
    double a, b, c, d;
}Four_doubles;

typedef struct marker{ // used in reading in data file
    int chromosome_number;
    int marker_number;
    int marker_position_on_chromosome; // 
    double marker_distance_on_chromosome; // distance from left end of chromosome
    double marker_distance; // distance from left end of genome
    char marker_name[MAXMARKERNAMELENGTH];
    char marker_sign[4]; // -, +, or ?
}Marker;

typedef struct marker_end* marker_end_ptr;

typedef struct marker_end{
    SHORTINT markend_numb; // 0 through 2*NMARKER+1    
    SHORTINT position; // position of the marker end in the permutation
    SHORTINT is_end_marker; // TRUE if marker_end is at end of chromosome (e.g. 0 and 2N+1 in 1 chromosome case)
    SHORTINT chromosome; // number of chromosome (0 - C-1) this marker end is on
    double distance; // distance from left end of whole genome
        // int arm; // LEFT or RIGHT, specifies whether marker is to L or R of centromere
    double d_black; // distance from adjacent marker ("length of black edge")
    double theta; // recombination fraction between this marker_end and neighbor
    int Nrecomb;  // the number of recombinations between this marker and neighbor, in some genotype data set
        //  char* name; // pointer to array of chars holding marker name. up to MAXMARKERNAMELENGTH-1 chars.
    marker_end_ptr black;
    marker_end_ptr other;
}Marker_end;

typedef struct permutation{
    Marker_end* p;      // array of Marker_end structures
    SHORTINT n_mark;    // number of real markers
    SHORTINT n_chrom; // number of chromosomes
    SHORTINT n_chrom_nonempty; //
    int n_inv, n_trans; // numbers of inversions, translocations ( n_inv + n_trans/2 == N+M choose 2 )
    double A_i, A_t; // Areas leading to inversions, translocations (A_i + A_t/2 == 1/2(L_g^2 - sum(l_i^2))
    SHORTINT n_markers_on_chromosome[MAX_N_CHROMOSOMES];
    double chromosome_lengths[MAX_N_CHROMOSOMES];
    int n_edge; // = n_mark + n_chrom
    double L_g; // genome length; (if no distances, use N+M as length)
    char** marker_names; // points to array of marker_names. Note these go from marker_names[0] to marker_names[n_mark-1]
}Permutation;

typedef struct reversal{
        // left and right are the positions of the marker ends at the
        // ends of the section to be reversed.
        // flip_chromosome is FALSE if right chromsome is not to be flipped
        // flip_chromosome is TRUE if right chromsome is to be flipped; but ignored if not a translocation
        // in that case the Reversal 's for the two translocations (flipped and not)
        // are the same except for the flip_chromosome field. the "right" field is
        // adjusted when the chromosome flip is done (in reverse_chromosome) so that
        // e.g. if  (1,2,3,4) and (5,6,7,8) are the two chromosomes, and
        // we flip the second one to get (8,7,6,5) and then translocate 3,4) and (8,
        // then left = 3, and right = 5, (the position of marker 8 after the flip)
        //  SHORTINT left, right;
    SHORTINT is_inversion, is_fission, same_eds;
    int numbs[4]; // numbs[i] are numbers of marker ends at ends of cut edge. All chromosome ends mapped to UNKNOWN
    int marker_numbers[4]; // like numbs, but without mapping of chromosome ends all to UNKNOWN
    int delta_c_i, delta_c_g; // delta c_i and delta c_g for this reversal
    int stuck; // if there was no delta_d = -1 inv/translocation available, this is TRUE
    double ld, rd, A; // distances of break points from L end of genome; A = l1*l2, i.e. prob of dist between a,b and c,d
    double dbreakod[4]; // distances of break points from marker_ends specified by marker_numbers, as fraction d_black
    int to_acbd; // TRUE -> marker_ends numbs[0],numbs[2] will be neighbors after rev, and also numbs[1],numbs[3]. FALSE -> numbs[0],numbs[3] neighbors, & numbs[1],numbs[2]
    int n_chroms_involved; // 1 for inversions, 2 for translocation.
    int chroms_involved[2]; // chromosome numbers of chromosome(s) involved in inv/trans
} Reversal;

typedef struct cycle_element{ // these correspond to black edges
    SHORTINT n1, n2;      // the numbers of the marker ends on either side of the black edge (in what order?)
        // SHORTINT p2chrom1, p2chrom2; // these are the chromosome numbers (in downstream perm, p2) of marker ends on either side of black edge 
    double d1, d2; // distances of marker_ends on either side of black edge from L end of genome
    double d12; // distance of marker_ends on either side from each other
    SHORTINT position;		// position = 0,(1, etc) for the black edge joining marker_ends in positions 0 and 1 (2 and 3)
    SHORTINT direction; // direction in which edge is crossed in going aroung the cycle, RIGHT or LEFT.
    int chromosome; // 0 to M-1, the number of the chromosome the element belongs to
    int eds; // eds number; numbered 0 to n_eds-1 where n_eds is the number of eds's in the whole cd
    int at_left_end; // TRUE if this is leftmost edge of chromosome
    int at_right_end; // TRUE if this is rightmost edge of chromosome
    struct cycle_element* next; // ptr to next cycle_element in cycle
} Cycle_element;

typedef struct cycle{ // these correspond to eds's i.e. "end delimited segments"
    SHORTINT n_black; // number of black edges in cycle
    SHORTINT n_gray;  // number of gray edges in cycle - just for checking in check_edss
    Cycle_element* start_elem_ptr; // ptr to beginning Cycle_element of cycle
    SHORTINT is_end_cycle;
    SHORTINT first_black, last_black; // first_black == TRUE <-> first edge on eds is a black edge.
} Cycle;

typedef struct cycle_decomposition{
    const Permutation* perm1, * perm2;
    SHORTINT cycle_found[2*NMARKERMAX + 2*MAX_N_CHROMOSOMES]; // TRUE if the cycle that marker end i is on has been found, else FALSE 
    Cycle_element elements[NMARKERMAX + MAX_N_CHROMOSOMES]; // array of cycle elements (black edges)
    Cycle edss[NMARKERMAX + MAX_N_CHROMOSOMES]; // array of eds's
    Cycle* end_edss[2*MAX_N_CHROMOSOMES]; // array of end eds's
    Cycle* int_edss[2*NMARKERMAX]; // array of int eds's
    int n_dci[3]; // number of reversals with delta c_internal = -1, 0, +1
    int n_dcg[3]; // number of reversals with delta c_g = -1, 0, +1
    int n_dd[3]; // number of reversals with delta d = -1, 0, 1. delta d = delta (c_g - c_i)
    int n_cycles; // total number of cycles in cycle decomposition
    int n_end_cycles; // total number of end cycles in cd
    int n_int_cycles; // number of internal cycles
    int n_mark;
    int n_chrom;
    int c_g; // number of eds's with gray edges on both ends (= number of eds's with black edges on both ends)
    int n_edss; // number of eds's in whole cycle decomposition
    int n_end_edss; // number of end eds's (doesn't count ones corresponding to empty chromosomes in target.
                    // n_end_edss should be equal to 2M-n_empty_chrom_in_target
    int n_empty_chrom_in_target; // number of empty chromosomes in target genome (perm2)
    int n00_seds, n00_deds1, n00_deds2; //
    int n00_seds_inv, n00_deds1_inv, n00_deds2_inv; //
    int n00_seds_trans, n00_deds1_trans, n00_deds2_trans; // 
    int d; // distance, = n_mark - n_chrom - n_int_cycles + c_g
    int n_dci_inv[3], n_dci_trans[3]; // number of inversions/translocations with delta c_internal = -1, 0, +1
    int n_dcg_inv[3], n_dcg_trans[3]; // number of inversions/translocations with delta c_g = -1, 0, +1
    int n_dd_trans[3], n_dd_inv[3]; // numbers of translocations, inversions with delta d = -1, 0, 1
    
    int n_dsw_inv_dci1[5], n_dsw_trans_dci1[5];
    
    int n_dsw_inv_dcgm1[5], n_dsw_trans_dcgm1[5];
    
    int n_dsw_dcgm1[5], n_dsw_dci1[5]; // sum of inv and trans
    int n_dsw_inv_ddm1[5], n_dsw_trans_ddm1[5];

    int n_ddm1_nsw2before[5];
    int n_ddm1_ncomdectozero[5];
    int n_ddm1_ncomdecrease[5];
    int n_ddm1_ncomnochange[5];
    int n_ddm1_ncomncorinc[5]; // no change or increase
    
} Cycle_decomposition;

typedef struct step{
    struct step* next, * prev; // pointer to next , previous, steps in a path
    Permutation* perm; // the permutation before the step is taken
    Reversal rev;  // the reversal which is applied to perm to get to the next permutation
    double prob; // for debugging
    int n_dum; // number of dummy events
}Step;

typedef struct path{
    Step* start;
    Step* end;
    double prob; // the probability of proposing this path from start to end (given the end points) (doesn't take account of dummies)
        // probably should eliminate this field (?) 
    int length; // number of inversions & translocations 
    int length_a; // number of steps on path (inversions, translocations and dummies)
        //  int length_r; // number of real events (inversions, translocations)
    int length_i; // number of inversions
    int length_t; // number of translocations on path
    int length_f; // number of fissions on path
    int L, L_i, L_t;
    int n_stuck; // number of times in generating path that no delta_d = -1 step is available
    int accept; // 1 if path resulted from accepted update, 0 otherwise
    int dd01_trans, dd01_inv; // TRUE if >0 delta_d = 0 or 1 step of spec. type (stuck cases not counted)
    double fission_factor; // = product over fissions of (2M_0) where M_0 is the number of empty chromosomes before fission
    double swap_end_factor;
}Path;

typedef struct lambdas{
    double lambdaI, lambdaT; // rates for inversions, translocations, respectively.
    double Lambda; // total rate, including dummy events (if any).
    double xi; //if LambdaI < lambdaT, xi = lambdaI/lambdaT  else xi = (2 - lambdaT/lambdaI)
    double r;
}Lambdas;

typedef struct state{
    Path* path;
    Lambdas* L;
    int Imax, Imin, Tmax, Tmin;
    int Mode; // 0 no uniformization, fixed Lambda, etc. 
}State;

typedef struct run_info_in{
    char* Genome_data_file; // name of file containing description of genomes
    int Use_distances; // Whether to make use of distance information (if present in the file)
    int N_rev_fake; // if > 0 is number of inversions/translocations to do to produce fake data; else use data from file as is
    double lambda_ratio_fake; // generate fake data with this value as ratio of inversion rate to translocation rate (lambdaI/lambdaT)
    int N_chain; // number of chains to run in parallel for data
  int N_short_chains; // number of chains to be initialized to short paths 
    int N_data; // number of simulated data to do
    
    int Update_lambdas_mode; // 0: no uniformization, fixed lambda 1: no uniformization, simulated Lambda
        // 2,3 lambdaI, lambdaT parametrization, 4,5,6: Lambda, r parametrization
    double Lambda; // For no-uniformization modes: <0 -> Lambda is simulated; if positive, is the fixed value of Lambda.
    double lambda_max; // prior p(lambdaI, lambdaT) = 1/lambda_max^2 for lambdaI, lambdaT both < lambda_max
    double Lambda_max; // Max value of Lambda. Value of Lambda corresponding to LambdaI = lambdaT = lambda_max
      
    double Epsilon; // Epsilon is parameter governing prob of picking delta c = +1, 0, -1 inversions, translocations.
    double Init_long_param, Init_long_sf_param; // these determine (approx) how long the long initial paths will be
        // epsilon_init is chosen such that long paths are approx init_long_param*sqrt(d) longer, and
        // sign flips are done with prob init_long_sf_param/sqrt(N), implying ~init_long_sf_param*sqrt(N) flipped signs (it this is <<N)

    int N_temperatures; // number of different temperatures for Metropolis-coupled MCMC (MCMCMC)
    double T_cold, T_hot; // Coldest and hottest temperatures. The N_temperatures are evenly spaced
    int Do_MC3_swap; // TRUE/FALSE, controls whether the swapping between different T chains is actually done (for checking purposes).
    int MC3_max_pathlength; // can use to limit path length for high T chains (which might blow up in length), by setting pi(path) to zero if L_path >this
  
//    int lambdaIT_update_mode; // 0: lambdaI_prop = lambdaI + lambdaIT_step_width_coeff*lambda_max
//    double lambdaIT_step_width_coeff; // 
//    double r_step_width; // r_prop = r + r_step_width*(drand() + 0.5)
    
    double Conv_srh; // value for sqrt(rhat) (using whole chains) at which conv is declared
    double Conv_srh_half2; // value for sqrt(rhat) (using 2nd half) at which conv is declared
    double Conv_Ldist; // find sigma of P(L|D) for each L (over chains), sum over L, divide by total number of steps*chains; when less than Conv_Ldist declare convergence
    double Conv_Lambdadist; // similar to Conv_Ldist, but using Lambda
    double Conv_rdist; // similar to Conv_Ldist, but using r
    double Conv_lambdaIdist; //
    double Conv_lambdaTdist; //
    int Burn_in; // number of steps of burn-in, or if < 0, burn-in is until convergence
    int Stop_at; // number of burn-ins to run. e.g. Stop_at = 5 -> if burn-in at N steps, run for 5N steps in all
    int Min_burn_in_length; // Go at least this many steps before declaring burn-in done
    int Max_runlength; // give up after this many iterations, even if no convergence.
     long Rand_seed; // RNG seed, controls which fake data generated
    int L_prop_dist_type; // controls which proposal distribution to use for lengths; 1: oldstyle, 2: 1-tanh, 3: lorentzian
    double Alpha; // location of peak (or falling edge of tanh) is at N*Alpha
    double Xi; // width of peak (or of transition region from hi to low) is alpha*N/Xi  
    int Choose_signs;  // 0 -> unsigned, 1,2,3,4, ... various algorithms for fiddling with signs for unsigned case.
    double Flip_one_sign_epsilon; // In range 0 to 1. If small, sign changes conducive to short paths are strongly preferred
  // so Epsilon, Flip_one_sign_epsilon, LambdaI_width_coeff, LambdaT_width_coeff, r_step_width,  Alpha, Xi are relevant to
  // 'size' of proposed step, and hence to acceptance prob.
  int Output_paths; // TRUE -> output file with paths
  double LambdaI_width_coeff; // determines size of proposed change in lambdaI
  double LambdaT_width_coeff; // determines size of proposed change in lambdaT
  double r_step_width; // determines size of proposed change in r (=lambdaI/lambdaT)

}Run_info_in;

typedef struct path_tree_node{
        //short revseq[MAXSHORTPATHLENGTH+1];
    short* revseq;
    Path* path;
    int count;
    double pi_rel, pi_rel_old;
    double fission_factor, swap_end_factor;
    struct path_tree_node* left; // pointer to left subtree containing items smaller than in this node
    struct path_tree_node* right; // pointer to left subtree containing items greater than in this node
}Path_tree_node;

typedef struct path_tree{
    Path_tree_node* start;
}Path_tree;

/* typedef struct revseq{ */
/*     short revarr[MAXSHORTPATHLENGTH + 1]; */
/* }Revseq; */

typedef struct stats{
    double N;
    double sum;
    double sumsqr;
    double mean;
    double stddev;
    double stddev_mean;
}Stats;

typedef struct histogram{
    int nbins;
    double minval; // lowest value which goes in bin zero
    double binwidth;
    double sum;
    double peak;
    int peak_bin;
    double mode; // center of peak bin
    int low_ne, hi_ne; // low_ne is smallest i s.t. histogram[i] > 0
    double* histogram;
    Stats* stats;
}Histogram;


typedef struct chain{
    State* state;
    double temperature;
    Histogram* L_hist, * r_inv_hist, * L_i_hist, * lambdaI_hist, * L_t_hist, * lambdaT_hist;
    int N_steps_so_far;
    int L_min, L_i_min, L_t_min; 
    double sumL, sumLsqr, W_L;
    double sumLi, sumLisqr, W_Li;
    double sumLt, sumLtsqr, W_Lt;
    double sumlambdaI, sumlambdaIsqr, W_lambdaI;
    double sumlambdaT, sumlambdaTsqr, W_lambdaT;
        // double n_prop_path, n_prop_lamI, n_prop_lamT, n_prop_r, n_prop_Lambda; 
    double n_acc_path, n_acc_lambdaI, n_acc_lambdaT, n_acc_r, n_acc_Lambda; 
}Chain;

typedef struct chain_set{
    int N_chain;
    int N_temperatures; 
    Chain*** chain; // array of N_chain pointers to arrays of N_temperature pointers to Chains
    Histogram* L_hist, * r_inv_hist, * L_i_hist, * lambdaI_hist, * L_t_hist, * lambdaT_hist;
    double Avg_L[MAX_N_TEMPERATURES], Avg_Li[MAX_N_TEMPERATURES], Avg_lambdaI[MAX_N_TEMPERATURES],
        Avg_Lt[MAX_N_TEMPERATURES], Avg_lambdaT[MAX_N_TEMPERATURES];
    double B_L[MAX_N_TEMPERATURES], B_Li[MAX_N_TEMPERATURES], B_lambdaI[MAX_N_TEMPERATURES],
        B_Lt[MAX_N_TEMPERATURES], B_lambdaT[MAX_N_TEMPERATURES]; 
    double W_L[MAX_N_TEMPERATURES], W_Li[MAX_N_TEMPERATURES], W_lambdaI[MAX_N_TEMPERATURES],
        W_Lt[MAX_N_TEMPERATURES], W_lambdaT[MAX_N_TEMPERATURES];
    double Pa_path[MAX_N_TEMPERATURES], Pa_lambdaI[MAX_N_TEMPERATURES], Pa_lambdaT[MAX_N_TEMPERATURES],
        Pa_r[MAX_N_TEMPERATURES], Pa_Lambda[MAX_N_TEMPERATURES];
}Chain_set;

// enum genotype {aa, aA, Aa, AA,  Aa_aA, ax_xa, Ax_xA, xx, ax, xa, Ax, xA, ChromEnd}; // In BC case, only AA, Aa, Ax occur
// enum genotype {AA, Aa, aA, aa,  Aa_aA, ax_xa, Ax_xA, xx, ax, xa, Ax, xA}; // In BC case, only AA, Aa, Ax occur




