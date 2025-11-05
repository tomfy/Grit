// histogram.h

Histogram* create_histogram(int nbins, double minval, double binwidth);
void free_histogram(Histogram* hist);
void insert_in_histogram(Histogram* hist, double value, double weight);
void print_histogram(FILE* fp, Histogram* hist, int normalize, double A);
void print_histograms(FILE* fp, int Nhists, Histogram** the_hists, int normalize, double A);

void initialize_stats(Stats* stat);
void insert_in_stats(Stats* stat, double value, double weight);
void print_stats(FILE* fp, Stats* stats);

void cred_int(Histogram* hist, double CL, double* xlo, double* xhi, double* ylo, double* yhi);
void cred_int1(double* x, double* a, int N, double CL, double* xlo, double* xhi, double* ylo, double* yhi);

