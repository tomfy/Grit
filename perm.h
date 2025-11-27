// perm.h

// Dealing with permutations

// read genomes from input file
// void get_data_from_file(FILE* fp, Permutation** genomes);
//int xxxy(FILE* fp, Permutation** genomes, int Use_distances);
//int get_data_from_file_new(FILE* fp, Permutation** genomes, int Use_distances);
int get_genomes_from_file(FILE* fp, Permutation** genomes, int Use_distances);
int get_genome_from_file(FILE* fp, Permutation** genome, int Use_distances);

Permutation* copy_perm(const Permutation* src); // return pointer to copy of permutation pointed to by src

// Dealing with reversals
int reverse(Permutation* perm, Reversal* rev);  // do specified reversal
Reversal make_signflip_rev(Permutation* perm, int cut1, int cut2); // given cut positions[0, N], return a reversal.
//Reversal choose_make_rev(int n_edges, Cycle_element* edge1, Cycle_element* edge2, int flip_chromosome);
Reversal make_rev(int n_edges, Cycle_element* edge1, Cycle_element* edge2, int flip_chromosome, const Permutation* perm);
Reversal make_rev1(const Permutation* perm, Marker_end* a, Marker_end* b, Marker_end* c, Marker_end* d, int flipped_chromosome);

void invert_translocate(Permutation* perm, Reversal* rev, int get_new_dists);
//int invert_translocate_old(Permutation* perm, Reversal* rev);
//int invert_translocate_new(Permutation* perm, Reversal* rev, int get_new_dists);

Reversal rand_reverse_r(Permutation* perm, double r, double* prob, Cycle_element* elements);
//int reverse_chromosome_old(Permutation* perm, Reversal* rev, int the_chromosome);
int reverse_chromosome(Permutation* perm, int the_chromosome, double* rd);
// Dealing with cycle decompositions
//inline
marker_end_ptr get_cycle_start(marker_end_ptr the_end, int* cycle_numbers);
int run_around_cycle(const Permutation* perm1, const Permutation* perm2, Cycle_decomposition* the_cd,
            marker_end_ptr cycle_start);
int get_CD(const Permutation* perm1, const Permutation* perm2, Cycle_decomposition* the_cd);
void get_cycle_lengths(Permutation* perm, Cycle_decomposition* the_cd);

int n_edge_pairs(Permutation* perm);
//int get_abcd(const Permutation* perm, const Reversal* rev, Marker_end** ap, Marker_end** bp, Marker_end** cp, Marker_end** dp);

int get_n_dcg1(Cycle_decomposition* the_cd);
int get_n_dcg1_old(Cycle_decomposition* the_cd);
int get_n_dcgm1(Cycle_decomposition* the_cd, const Permutation* p2);
int get_n_dci1(Cycle_decomposition* the_cd, const Permutation* p2);
int get_n_dci1_eds(int N_chrom, Cycle_element** start, int* n_dci0, int* n_dci1_inv, int* n_dci1_trans);
int get_n_dcim1(Cycle_decomposition* the_cd);
void get_end_int_edss(Cycle_decomposition* the_cd);
int get_n00_disteds2(Cycle_decomposition* the_cd);
int get_c_g(Cycle_decomposition* the_cd);

int count_between_chrom_breaks(const Cycle_element* edge1, const Cycle_element* edge2, const Permutation* p2);
int delta_n_switches(const Cycle_element* edge1, const Cycle_element* edge2, const Permutation* p2, int* flipped_result,
                       int* nswbefore, int* nswafter, int* nswafterflipped);
int delta_n_switches1(const int* marker_numbs, const Permutation* p1, const Permutation* p2, int* flipped_result);

Reversal get_spec_dci1_rev(Cycle_decomposition* the_cd, int spec_rev, int type);
Reversal get_spec_dci1_rev_new(Cycle_decomposition* the_cd, int spec_rev, int type, int spec_dsw_index);
Reversal get_spec_dci1_rev_new1(Cycle_decomposition* the_cd, int spec_n, int type, int spec_dnsw_index);
Reversal get_spec_dcim1_rev(Cycle_decomposition* the_cd, int n, int type);
Reversal get_spec_dcgm1_rev(Cycle_decomposition* the_cd, int n, int type);
Reversal get_spec_dcgm1_rev_new(Cycle_decomposition* the_cd, int spec_n, int type, int spec_dsw_index);
Reversal get_spec_dcg1_rev(Cycle_decomposition* the_cd, int n, int type);

Reversal get_spec_00seds_rev(Cycle_decomposition* the_cd, int n_spec, int type);
Reversal get_spec_00deds1_rev(Cycle_decomposition* the_cd, int n_spec, int type);
Reversal get_spec_00deds2_rev(Cycle_decomposition* the_cd, int n_spec, int type);


// void get_data_from_file_new(FILE* fp, Permutation** genomes);
int
marker_cmp(const void* m1, const void* m2);
void
print_perm_array(int Nchrom, int* Nmarkers, Marker** perm_array);
int
same_markers(int N, Marker* g1, Marker* g2);


double
get_delta_d_prob_new(Cycle_decomposition* the_cd, int delta_d, int same_eds, int is_inversion, int n_inv, int n_trans, double Epsilon, double lambda_ratio,
                      int dnsw);
Reversal
delta_ci_cg_reverse_new(Cycle_decomposition* the_cd, double epsilon, int n_inv, int n_trans, double* prob, double lambda_ratio);
double
p_dnsw(int dnsw_index, int n);

void get_chroms_for_switches(int n0, int n1, int ab, int* chrom_a, int* chrom_b, const Permutation* p1,  const Permutation* p2);
int delta_n_switches2(const int* marker_numbers, const Permutation* p1, const Permutation* p2, int* flipped_result);
int n_switches(const int* marker_numbers, const Permutation* perm, const Permutation* targ_perm);
int n_switches1(const int* marker_numbers, const Permutation* perm, const Permutation* targ_perm);
int n_syn(const Permutation* p1, const Permutation * p2);
int chroms_in_common(Cycle_element* edge1, Cycle_element* edge2, const Permutation* p1, const Permutation* p2,
                     int* ncom_acbd, int* ncom_adbc);
//int get_abcd(const Permutation* perm, const Reversal rev, Marker_end** ap, Marker_end** bp, Marker_end** cp, Marker_end** dp);
//int get_abcd_old(const Permutation* perm, const Reversal rev, Marker_end** ap, Marker_end** bp, Marker_end** cp, Marker_end** dp); 
//int get_abcd_new(const Permutation* perm, const Reversal* rev, Marker_end** ap, Marker_end** bp, Marker_end** cp, Marker_end** dp);
//Marker_end* get_me_pair(const Permutation* perm, int* pos, Marker_end** e, Marker_end* prev_e_empty);

int get_4ends(const Permutation* perm, int* numbs, Marker_end** ends, Marker_end** ordered_ends, int* order);
Marker_end* fix_me_pair(const Permutation* perm, Marker_end** e, Marker_end* prev_e_empty);
int order_4ends(Marker_end** ends, Marker_end** ordered_ends, int* order);
void get_signflip_break_dists(Marker_end* a, Marker_end* b, Marker_end* c, Marker_end* d, double* dbreakod, double* ld, double* rd);
void get_rand_break_dists(Marker_end* a, Marker_end* b, Marker_end* c, Marker_end* d, double* dbreakod, double* ld, double* rd);

int check_rev(int n_markers, Reversal* rev);

int invert_translocate_new1(Permutation* perm, Reversal* rev, int get_new_dists);
int invert_translocate_new2(Permutation* perm, Reversal* rev, int get_new_dists);
int reverse_new2(Permutation* perm, Reversal* rev);

void
fix_rev1(Permutation* perm, Reversal* rev);

Marker_end* get_chromosome_left_end(Marker_end* m);

void simple_rand_reverse(Permutation* perm, double r);
void print_oxford_grid(FILE* fp, Permutation* perm1, Permutation* perm2);

// end of perm.h

