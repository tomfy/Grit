// path_etc.h   declarations of functions for allocation, checking, printing, etc.
// related to paths

// Memory allocation and deallocation, checking, printing
void freepath(Path* path);
Path* path_alloc(Step* step);
Step* step_alloc(Permutation* perm);

// checking
int check_path_consistency_multi(const Path* the_path);
int check_path_lengths(const Path* path);

double test_hits_vs_probs_1_step(FILE* fptest, Permutation* p1, Permutation* p2, double Epsilon, int n_tries, int check_every);
//double test_path_hits_vs_probs(FILE* fptest, Permutation* p1, Permutation* p2, double Epsilon, int n_tries, int check_every);

// print info
void print_dummies(Path* path);
void print_rev(FILE* fp, Reversal rev);
void print_rev_brief(FILE* fp, Reversal rev);
void print_path(FILE* fp, const Path* the_path);
void print_path_brief(FILE* fp, const Path* the_path);
//void print_path_trans(FILE* fp, Path* the_path);

double get_fission_factor(const Path* path);
double get_swap_end_factor(const Path* path);

double f_Ltheat(double T_hot, double T_inv, int Lt);
double f_Ltheat1(double T, double beta, double Lambda, double lambdaT, int Lt);
void teststates(Permutation* p);

void print_dummies1(FILE* fp, const Path* the_path);
void print_dummies2(Step* step);

int check_path_revs_ordered(Path* path);


int check_rev_perm_consistent(Permutation* perm, Reversal* rev);
int rev_just_swaps_chromosome_ends(Permutation* perm, Reversal rev);
// end of path_etc.h
