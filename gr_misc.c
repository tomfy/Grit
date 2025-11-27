// gr_misc.c   miscellaneous functions used in grit

#include "grit.h"
#include "gr_misc.h"
#include "float.h"
extern int Error_count;

// sort an array
void qqsort(int v[], int left, int right)
{
	int i, last;
	if ((right - left) <= 10) { // array with 10 or fewer elems are sorted by straight insertion
		insertionsort(v+left, 1+right-left); // use n^2 method if array issmall
		return;
	}
	if (left >= right) // v has only 1 element - no sorting needed
		return;

	// is the following line needed?
	// Yes, avoids n^2 behavior in worst (pre-ordered) case.
  swap(v, left, (left+right)/2);
	last = left;
	for (i = left+1; i <= right; i++)
		if(v[i] < v[left])
			swap(v, ++last, i);
	swap(v, left, last);
	qqsort(v, left, last-1);
	qqsort(v, last+1, right);
}
void swap(int v[], int i, int j)
{
	int temp;
	temp = v[i]; v[i] = v[j]; v[j] = temp;
}

void insertionsort(int v[], int size)
{
	int i, j, temp;
	for (j=1; j<size; j++){
		temp = v[j];
		for (i=j-1; i>=0; i--){
			if (v[i] <= temp) goto abc;
			v[i+1] = v[i];
		}
		abc: v[i+1] = temp;
	}
}

double factorial(int a)
{
    if(a < 0){
        printf("arg of factorial function is negative: %i \n", a);
        return -1.0;
    }
    if(a <= 1) return 1.0;
    else {
        return (double)a*factorial(a-1);
    }
}

double log_factorial(int a)
{
        // calculates log(a!) recursively as log(a) + log((a-1)!)
    if(a < 0){
        printf("arg of factorial function is negative: %i \n", a);
        return -1.0;
    }
    if(a <= 1) return 0.0;
    else {
        return log((double)a) + log_factorial(a-1);
    }
}

double step_fcn(int i) 
{ 
        // returns 1.0 if i>= 1, 0.0 otherwise; 
    if(i >= 1) return 1.0; 
    else return 0.0; 
}  // end of function step_fcn

int
compare_int(int* a1, int* a2)
{
    if(*a1 > *a2) return 1;
    else if(*a2 > *a1) return -1;
    else return 0;
}


int
int_min(int i, int j)
{
    int temp;
    temp = i<j? i:j;
    return temp;
}


double frac_diff(const double a, const double b)
{
    return (a - b)/(a + b);
}

int
real_close(double a, double b, double threshold)
{
        // returns 1 (TRUE) if fractional difference of args is < threshold
        // 0 otherwise
    double x;

        //   printf("in real_close, threshold: %g \n", threshold);
    if(threshold > 1.0e-7) { printf("in real_close. threshold > 1.0e-7. threshold : %g \n", threshold); getchar(); }
    if(a == 0.0 && b == 0.0) return 1;
    x = fabs(frac_diff(a, b));
    if(x < (threshold)){ return 1; }
    else {
        printf("In real_close; fractional difference of arguments is too big. \n arguments: %22g %22g  fabs((a-b)/(a+b)): %g, threshold: %g \n",
               a, b, x, threshold);
        return 0;
    }
} // end of real_close

int
real_close2(double a, double b, double threshold, double Q)
{
    
    double x;

        //   printf("in real_close, threshold: %g \n", threshold);
    if(threshold > 1.0e-7) { printf("in real_close2. threshold > 1.0e-7. threshold : %g \n", threshold); getchar(); }
    if(a == 0.0 && b == 0.0) return 1;
    x = fabs(a-b)/(Q + fabs(a+b));
    if(x < (threshold)){ return 1; }
    else {
        printf("In real_close2; arguments not close enough. \n arguments: %22g %22g  |a-b|/(Q + |a+b|): %g, threshold: %g , Q: %g\n",
               a, b, x, threshold, Q);
        return 0;
    }
} // end of real_close2

int
logs_real_close(double loga, double logb, double threshold)
{
    double x, denom;

    if(threshold > 1.0e-7) { printf("in logs_real_close. threshold : %g \n", threshold); getchar(); }
    x = fabs(loga - logb);
    denom = 0.5*fabs(loga + logb);
    if(denom < 1.0) denom = 1.0;
        //  if(x/denom > 1.0e-14) printf("in logs_real_close: %g %g %g \n", loga, logb, x/denom);
    if(x/denom < (threshold)){ return 1; }
    else {
        printf("In logs_real_close; difference of arguments is too big. \n loga, logb: %22g %22g  loga-logb %g, threshold: %g \n",
               loga, logb, x, threshold);
        return 0;
    }
} // end of real_close 


double** 
alloc_double2d(int n, int m)
{
        // allocates
        // 1) an array of n pointers to arrays of m doubles,
        // 2) the n arrays of m doubles.
        //   initialized to zero?
    double** x;
    int i;

    x = (double**)chcalloc((size_t)n, sizeof(double*));
    for (i=0; i<n; i++){
        x[i] = (double*)chcalloc((size_t)m, sizeof(double));
    }
    return x;    
} // end of alloc_double2d

int** 
alloc_int2d(int n, int m)
{
        // allocates
        // 1) an array of n pointers to arrays of m doubles,
        // 2) the n arrays of m ints.
        //   initialized to zero?
    int** x;
    int i;

    x = (int**)chcalloc((size_t)n, sizeof(int*));
    for (i=0; i<n; i++){
        x[i] = (int*)chcalloc((size_t)m, sizeof(int));
    }
    return x;    
} // end of alloc_int2d



void 
free_double2d(int n, double** x)
{
        // free an array of n pointers to arrays of doubles
    int i;
    for (i=0; i<n; i++){
        free(x[i]);
    }
    free(x);
} // end of free_double2d

void 
free_int2d(int n, int** x)
{
        // free an array of n pointers to arrays of ints
    int i;
    for (i=0; i<n; i++){
        free(x[i]);
    }
    free(x);
} // end of free_int2d

int
is_normal(double x)
{
        // returns TRUE iff x is "normal" , i.e. not inf, NaN, denormal or zero
        //  there is supposed to be a built in isnormal macro, and an FP_NORMAL macro
        // but I get errors during compilation when I use either of them.
        // return (__fpclassify(x) == MYFPNORMAL);

        // the above didn't work with M$ compiler, just ask that number be >0 and not too small or big
    return (x > 1.0e-300 && x < 1.0e300);
}

int
is_normal_or_zero(double x)
{
        // returns TRUE iff x is "normal" or zero, i.e. not inf, NaN, denormal
        //  there is supposed to be a built in isnormal macro, and an FP_NORMAL macro
        // but I get errors during compilation when I use either of them.
        //   return ((__fpclassify(x) == MYFPNORMAL) || (x == 0.0));

        // above didn't work with MS compiler.
        // I use this to test logs; just test they are not too small or too big
        // i.e. primarily not +-inf, or nan. some subnormal (close to zero) will
        // give true, others (close to inf) will give false
    return(x > -1.0e300 && x < 1.0e300);
}

void
check_for_null_pointer(void* p)
{
    if(p == NULL){
        printf("malloc returned NULL pointer. \n");
        process_error();
    }
}

void* chmalloc(size_t size)
{ //malloc with very simple checking for NULL pointer
    void* p = malloc(size);
    if(p != NULL) return p;
    printf("malloc returned NULL pointer. \n");
    process_error();
}

void* chcalloc(size_t n, size_t element_size)
{  // calloc with very simple checking for NULL pointer
    void* p = calloc(n, element_size);
    if(p != NULL) return p;
    printf("calloc returned NULL pointer. \n");
    process_error();
}

/* void double_swap(double* x, double* y){double temp = *x; *x = *y; *y = temp;} */
/* void int_swap(int* x, int* y){int temp = *x; *x = *y; *y = temp;} */
/* double beta(double theta, int R, int N){return pow(theta, R)*pow(1.0 - theta, (N-R)); } */

void process_error(void)
{
        // call this whenever an error is encountered
        // counts errors and exits when too many
    Error_count++;
    printf("Error_count: %i \n", Error_count);
    if(Error_count >= ERROR_COUNT_LIMIT){
        fprintf(stdout, "Exitting after %i errors.\n", Error_count);
        fflush(NULL);
        exit(EXIT_FAILURE);
    }
}

double beta_inc_half_old(int an, int aN)
{
        // integral from 0 to 1/2 of x^an*(1-x)^(aN-an)dx
        // this is special case of incomplete beta function
    static int N = -1;
    static double* resarray = NULL;
    static double* resarray2 = NULL;

    if(N != aN || N < 0){ //
        if(N < 0){ /* printf("first time calling beta_inc_half. n, m, N: %i %i %i\n", an , aN - an , aN ); */}
        else { printf("calling beta_inc_half; N changed. old N, new N: %i %i\n", N, aN );}
        
        int n, m;
        double dN = (double)aN;
        N = aN;
        if(resarray != NULL) free(resarray);
        resarray = (double*)chcalloc((N+1), sizeof(double));
        
        resarray[0] = 1.0/(dN + 1.0); // result for x^N*(1-x)^0 *2^N+1
        for(m=1, n = N-1; m<=N; m++, n--){  
            resarray[m] = (1.0 + (double)m*resarray[m-1])/(double)(n+1);
        }

        double denom;
        if(resarray2 != NULL) free(resarray2);
        resarray2 = (double*)chcalloc((N+1), sizeof(double));
        denom = pow(2.0, (double)(N+1));
        for(m=0; m<=N; m++){  
            resarray2[m] = resarray[m]/denom;
        }
            //  printf("bottom of beta_inc_half first time. \n");
    }
        //  printf("beta_inc_half. n, m, N: %i %i %i \n", an, aN, aN - an);
        //   printf("beta_inc_half return value: %g \n", resarray2[aN - an]);
        //  return resarray[aN - an]/pow(2.0, (double)(aN+1));
    return resarray2[aN - an];
} // end of beta_inc_half


double beta_inc_half(int n, int N)
{
        // integral from 0 to 1/2 of x^n*(1-x)^(N-n)dx
        // this is special case of incomplete beta function
    static int Nmax;
    static int first = TRUE;
    static double** resarray = NULL;
    static double** resarray2 = NULL;
    int i;
  
    if(first){ //
        Nmax = N;
        resarray = (double**)chcalloc((Nmax+1), sizeof(double*));
        resarray2 = (double**)chcalloc((Nmax+1), sizeof(double*));
        for(i=0; i<= Nmax; i++){
            resarray[i] = resarray2[i] = NULL;
        }
            // printf("first time calling beta_inc_half. n, m, N: %i %i %i\n", n , N - n , N );
        first = FALSE;
    }
    else if(N > Nmax){
        printf("In beta_inc_half. N, Nmax: %i %i \n", N, Nmax);
        printf("arguments to beta_inc_half. n, N: %i %i \n", n, N);
        getchar();
    }
    
   
    if(resarray[N] == NULL){ //
        int n, m;
        double denom;
        resarray[N] = (double*)chcalloc((N+1), sizeof(double));
        resarray2[N] = (double*)chcalloc((N+1), sizeof(double));
        denom = pow(2.0, (double)(N+1));
        resarray[N][0] = 1.0/((double)N + 1.0); // result for x^N*(1-x)^0 *2^N+1
        resarray2[N][0] = resarray[N][0]/denom;
        for(m=1, n = N-1; m<=N; m++, n--){  
            resarray[N][m] = (1.0 + (double)m*resarray[N][m-1])/(double)(n+1);
            resarray2[N][m] = resarray[N][m]/denom;
        }
    }
       
    return resarray2[N][N - n];
} // end of beta_inc_half



double beta_density_0half(int R, int N, double theta)
{
        // if theta is distributed with density
        // proportional to theta^R * (1-theta)^(N-R)
        // between 0 and 0.5
        // returns normalized prob density at theta.
        // so integral of beta_density_0half(R, N, theta) from 0 to 0.5 should be 1
    return pow(theta, (double)R)*pow(1.0 - theta, (double)(N - R))/beta_inc_half(R, N);
} // end of beta_density_0half

double theta_to_haldane(double theta, double inv_theta_knee)
{
    return -0.5*log(1.0 - 2.0*(theta/(1.0 + theta*inv_theta_knee)));
}


//inline
void double_swap(double* x, double* y){double temp = *x; *x = *y; *y = temp;}
//inline
void int_swap(int* x, int* y){int temp = *x; *x = *y; *y = temp;}
//inline
double beta(double theta, int R, int N){return pow(theta, R)*pow(1.0 - theta, (N-R)); }

// end of gr_misc.c

