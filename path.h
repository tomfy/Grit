// path.h

// Generating path from one permutation to another
void append_step_to_path(Path* path, Step* step);
Reversal take_step_multi(Permutation* perm, Cycle_decomposition* the_cd, double Epsilon, double* step_prob, double lambda_ratio);
Path* gen_path(Permutation* p1, Permutation* p2, double Epsilon, double lambda_ratio, double limit_l, double* log_prob); // generate a path from p1 to p2
//Path* gen_new_dist_path(const Step* olddownstart, const Step* newsubpathend);

void add_dummies(Path* path, const Lambdas* L); // add dummy events for uniformization
//inline
double lambdum_over_Lambda(const Permutation* p, const Lambdas* L);
//inline
double Lambda_real_over_Lambda(const Permutation* p, const Lambdas* L);

// Updating path
void get_section(Path* the_path, int update_length, /*  Step** before_step,  */Path* subpath);
double get_length_prob(const Run_info_in* r_in, int N, int L, int length);
void get_rand_length(const Run_info_in* r_in, int N, int L, int* length, double* prob);

// double prod_a_to_n(const Path* path, const Lambdas* L);
double log_prod_a_to_n(const Path* path, const Lambdas* L);
double log_prod_A(const Path* path);
double prod_1minusa(const Path* path, const Lambdas* L);
double get_n_inv_avg(const Path* path);
double get_n_trans_avg(const Path* path);

//double pi_ratio(int l_old, int l_t_old, int l_new, int l_t_new, int NpM, double lambda);  // used by the lambdaT = 0.5*lambdaI algorithm
//double pi_ratio_1(int l_old, int l_new, double lambda); // does just the exp(-Lambda)*Lambda^L/L! part of ratio
//double pi_ratio_check(int l_old, int l_t_old, int l_new, int l_t_new, int NpM, double lambda);
//double pi_ratio_indep(const Path* old_path, const Path* new_path, const int NpM, const Lambdas* L);  // used independent lambdas case

double log_pi_ratio(int l_old, int l_t_old, int l_new, int l_t_new, int NpM, double Lambda);  // used by the lambdaT = 0.5*lambdaI algorithm
double log_pi_ratio_1(int l_old, int l_new, double lambda); // does just the exp(-Lambda)*Lambda^L/L! part of ratio
double log_pi_ratio_check(int l_old, int l_t_old, int l_new, int l_t_new, int NpM, double Lambda);
double log_pi_ratio_indep(const Path* old_path, const Path* new_path, const int NpM, const Lambdas* L);

int get_path_lengths(Path* the_path);

int update_path(State* state, const Run_info_in* r_in, double exponent, int update_to_end); //
//int update_path_o(State* state, const Run_info_in* r_in, double exponent, int update_to_end); //
int update_path_2ways(State* state, const Run_info_in* r_in, double exponent);
int compare_revs(Reversal rev1, Reversal rev2);
//int rev_just_swaps_chromosome_ends(Reversal rev);
//int rev_just_swaps_chromosome_ends(Permutation* perm, Reversal rev);
double get_epsilon(Cycle_decomposition* the_cd, double Epsilon, int path_length, double limit_l);
double get_path_prob_multi(const Path* the_path, const Permutation* targ_perm, const double Epsilon, double lambda_ratio, double* log_prob) ;

int compare_paths(Path* path1, Path* path2, int* count);
Path* copy_path(Path* src_path);

double get_step_prob_multi(Permutation* p, Permutation* next_p, Permutation* targ_p, Reversal* rev, double Epsilon, double lambda_ratio);
double get_sum_probs(Permutation* p1, Permutation* p2, double Epsilon);

double pi1(const State* state);

double rel_pi(State* state);

Path* shorten_path(Path* path);
int find_whether_same_eds(Permutation* perm, Cycle_decomposition* the_cd, Reversal* rev);

double pi_ratio_gen(Path* p_old, Lambdas* L_old, Path* p_new, Lambdas* L_new);
double log_pi_ratio_gen(Path* p_old, Lambdas* L_old, Path* p_new, Lambdas* L_new);
double pi_ratio_gen1(int old_length, int old_length_a, int old_length_i, int old_length_t, double old_Lambda, double old_lambdaI, double old_lambdaT,  
                         int new_length, int new_length_a, int new_length_i, int new_length_t, double new_Lambda, double new_lambdaI, double new_lambdaT);
double log_pi_ratio_gen1(int old_length, int old_length_a, int old_length_i, int old_length_t, double old_Lambda, double old_lambdaI, double old_lambdaT,  
                         int new_length, int new_length_a, int new_length_i, int new_length_t, double new_Lambda, double new_lambdaI, double new_lambdaT);  

void get_LLiLt(Path* path);
int check_LiLt(const Path* path);
void do_zeta_hist(const Path* path);

void track_lengths(const Path* path, double* the_ptl, double* the_atl);
void possible_track_lengths_1chrom(int chrom_size, double* the_ptl);

double get_A_i_avg(const Path* path);
double get_A_t_avg(const Path* path);



// int check_rev_perm_consistent(Permutation* perm, Reversal* rev);

void
continue_path_with_new_distances(Path* newpath, Step* oldsubpathend, double* logprodAold, double* logprodAnew);

Path* regenerate_path(Path* old_path);
int check_path_selfconsistent(Path* the_path);
void print_path_chroms_involved(const Path* path);


void update_path_order(Path* path);
void commutation_path(Path* path, Step* step, Step** start, Step** end, int* backward, int* forward);
int revs_commute(Reversal* rev1, Reversal* rev2);
int
check_path_next_prev(Path* path);
void
print_path_next_prev(Path* path);
// end of path.h
