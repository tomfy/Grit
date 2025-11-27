// perm_etc.h

int get_n_inv(Permutation* perm); 
int get_n_inv1(Permutation* perm);
int get_n_trans(Permutation* perm);
int get_n_inv_max(Permutation* perm);
int get_n_inv_min(Permutation* perm);
int get_n_trans_max(Permutation* perm);
int get_n_trans_min(Permutation* perm);

void
set_Ai_At(Permutation* perm);

// Checking functions
int compare_interior_edges(Permutation* perm1, Permutation* perm2, int* N1, int* N2);
int are_perms_equivalent(const Permutation* perm1, const Permutation* perm2);
int are_perms_equivalent_fast(Cycle_decomposition* the_cd, Permutation* perm1, Permutation* perm2);
int are_perms_equivalent_both(Cycle_decomposition* the_cd, Permutation* perm1, Permutation* perm2);
int are_perms_equal(const Permutation* p1, const Permutation* p2);
int check_markend_numbs(const Permutation* perm); // checks whether markend_numb fields match array index
int check_markend_is_end(const Permutation* perm);
int check_n_chrom_nonempty(const Permutation* perm);
int check_markend_positions(const Permutation* perm);
int check_directions(const Permutation* perm, Cycle_decomposition* the_cd);
int check_markend_black_ptrs(const Permutation* perm);
int check_markend_other_ptrs(const Permutation* perm);
int check_perm_selfconsistent(const Permutation* perm);
int check_markend_dists(const Permutation* perm);
int check_markend_dists2(const Permutation* perm);
int check_Ai_At(Permutation* perm);
void check_edss(Cycle_decomposition* the_cd);

// allocation, initialization
void permfree(Permutation** perm);
Permutation* init_genome_to_random(int N_chromosomes, int N_markers, int N_randomize);
Permutation* multipermalloc(int n_chromosomes, int n_markers);  // allocates memory for Permutation - multiple chromosomes version
Permutation* init_genome(int N_chromosomes, int* N_markers, double* chromosome_length, Marker perm_array[][MAX_N_MARKERS_PER_CHROMOSOME], double L_g);
void init_to_given_multi_perm(int* N_markers, double* chromosome_length, Marker in[][MAX_N_MARKERS_PER_CHROMOSOME], Permutation* out);
//void initialize_recombination_fractions(Permutation* perm);

// print info
void print_markers(const Permutation* perm);
void print_edges(Cycle_decomposition* the_cd, const int n_edges);
void print_perm(FILE* fp, const Permutation* p);
void print_perm_brief(FILE* fp, const Permutation* p);
void print_perm1(FILE* fp, const Permutation* p, int Nindiv);
void print_perm2(FILE* fp, const Permutation* p);
char* marker_end_name(const Permutation* p, Marker_end* m);
void print_perm_markers(FILE* fp, const Permutation* perm);
void print_mends(const Permutation* perm);
//void print_cycles(Cycle_decomposition* the_cd);
void print_edss(Cycle_decomposition* the_cd);

int check_rev_ordered(const Permutation* perm, Reversal rev);



/// end of perm_etc.h
