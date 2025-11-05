// chain?.h

struct XBW{
  double X, B, W, X2, B2, W2;
};

void check_run_params(const Run_info_in* r_in);
int converged(int* lengths, int fresh_chains);
int X_converged(int index, const Run_info_in* r_in,  double** cumeSteps, double** cumeX, double** cumeXsqr, struct XBW* the_XBW);
//void update_lambda(State* state, const Run_info_in* r_in);
Four_doubles update_lambdas(State* state, const Run_info_in* r_in, double exponent);
inline int update_lambdaI(const Path* path, Lambdas* L, double width, double lambda_max, double exponent);
inline int update_lambdaT(const Path* path, Lambdas* L, double width, double lambda_max, double exponent);
inline int update_r(Path* path, Lambdas* L);
inline int update_xi(Path* path, Lambdas* L);
void update_chain(Chain* the_chain, const Run_info_in* r_in);
void initialize_lambdas(State* state, const Run_info_in* r_in);
//double get_Lambda1(int Imin, int Imax, int Tmin, int Tmax, double lambdaI, double lambdaT);
double get_Lambda_old(int n_edges, double lambdaI, double lambdaT);
//void get_lambdaIT_from_Lambda_r_alt(int n_edges, double Lambda, double r, double* lambdaI, double* lambdaT);

double get_xi(double lambdaI, double lambdaT);
void get_lambdaIT_from_Lambda_xi_old(int n_edges, double Lambda, double xi, double* lambdaI, double* lambdaT);
void get_lambdaIT_from_Lambda_r_old(int n_edges, double Lambda, double r, double* lambdaI, double* lambdaT);

double get_Lambda(double L_g, double lambdaI, double lambdaT);
void get_lambdaIT_from_Lambda_xi(double L_g, double Lambda, double xi, double* lambdaI, double* lambdaT);
void get_lambdaIT_from_Lambda_r(double L_g, double Lambda, double r, double* lambdaI, double* lambdaT);
//int run_chains(Permutation* p1, Permutation* p2, const Run_info_in* r_in, double* cpu_time, int Chain_mode);

State** get_init_states(Permutation* p1, Permutation* p2, const Run_info_in* r_in, int n_short, int n_long, int max_path_length);

Chain_set* construct_chain_set(const Permutation* p1, const Permutation* p2, const Run_info_in* r_in);
Chain* construct_chain(const Permutation* p1in, const Permutation* p2in, const Run_info_in* r_in, int make_short_path, double Temperature);


void output_run_params(FILE* fp, const Run_info_in* r_in);
void input_run_params(FILE* fp, Run_info_in* r_in);

// path tree stuff went to path_tree_etc.h

void output1(int N_chain, State** the_state);
void output2(int nsteps, struct XBW ABW, struct XBW LBW, struct XBW lambdaBW);
void output4(FILE* fp, Path_tree* path_tree, Histogram* L_hist, Histogram* Lambda_hist,
             Histogram* L_i_hist, Histogram* lambdaI_hist, Histogram* L_t_hist, Histogram* lambdaT_hist);
void output5(FILE* fp, int steps_so_far, int i, double** cumeL, double** cumeLambda, int* Lengths, int N_chain);
void output5_old(FILE* fp, int steps_so_far, int i, double** cumeL, int N_chain);
void output7(FILE* fp, int Nchains, Histogram** L_hist, Histogram** lambda_hist);
void print_Lambdas(Lambdas* L);
double dist_converge(int N_chain, Histogram** hist);

int path_array_insert(Path* path_array[], int count_array[], Path* path); // like path tree but do it dumb way with array

double test_trans_sym(int size, long trans[][100]);

double srh(double a);

double check_path_freq_vs_pi_old(Path_tree* tree);

inline double within_chain_variance(int N, double sumX, double sumXsqr);
void chain_set_abw(Chain_set* chain_set);

// end of chain.h


