// run_chains.h

long run_chains(Chain_set* the_chain_set, const Run_info_in* r_in, long* cpu_time);

int fixed_lambda_mode(int mode);
Histogram* make_scaled_histogram(double xmin, double xrange, double minbw, double* binwidth);
int get_histogram_parameters(double xmin, double nbins, int Integers, double* xmin_out, double* bw);

double BoW(double B, double W);
long progress_output(const Run_info_in* r_in, int* n_mc3_swap, Chain_set* the_chain_set,
		     long steps_so_far, long n_mc3_swap_try, long tstart);
