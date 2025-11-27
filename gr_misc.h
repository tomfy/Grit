// gr_misc.h declarations of various general purpose functions

void qqsort(int v[], int left, int right);
void swap(int v[], int i, int j);
void insertionsort(int v[], int size);

double factorial(int a);
double log_factorial(int a);
int compare_int(int* a1, int* a2);
double step_fcn(int i);

int int_min(int i, int j);
double frac_diff(const double a, const double b);

int real_close(double a, double b, double threshold); //test whether a and b are close to equal (looks at abs((a-b)/(a+b)) )
int real_close2(double a, double b, double threshold, double Q);
int logs_real_close(double loga, double logb, double threshold);

double** alloc_double2d(int n, int m);
void free_double2d(int n, double** x);

int** alloc_int2d(int n, int m);
void free_int2d(int n, int** x);

int is_normal(double x);
int is_normal_or_zero(double x);

void check_for_null_pointer(void* p);
void* chmalloc(size_t size);
void* chcalloc(size_t n, size_t element_size);

void double_swap(double* x, double* y); // {double temp = *x; *x = *y; *y = temp;}
void int_swap(int* x, int* y); // {int temp = *x; *x = *y; *y = temp;}
double beta(double theta, int R, int N); //{return pow(theta, R)*pow(1.0 - theta, (N-R)); }

void process_error(void);

double beta_inc_half_old(int an, int aN);
double beta_inc_half(int an, int N); 
double beta_density_0half(int R, int N, double theta);
double theta_to_haldane(double theta, double inv_theta_knee);
