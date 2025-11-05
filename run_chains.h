// run_chains.h


int run_chains(Permutation* p1, Permutation* p2, const Run_info_in* r_in, double* cpu_time);
int run_chains_new(Chain_set* the_chain_set, const Run_info_in* r_in, double* cpu_time);

int fixed_lambda_mode(int mode);
Histogram* make_scaled_histogram(double xmin, double xrange, double minbw, double* binwidth);
int get_histogram_parameters(double xmin, double nbins, int Integers, double* xmin_out, double* bw);

double BoW(double B, double W);
