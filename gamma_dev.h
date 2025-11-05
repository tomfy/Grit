// gamma_dev.h
// 

double loggamma(double xx);
double gseries(double a, double x);
double gcf(double a, double x);
double incomplete_gamma(double a, double x);
//double gamma_dev3(double a);
double gamma_dev2(double a, double xmax);
double gamma_dev1(int k, double xmax); // based on routine in Numerical Recipes 
double gamma_dev(int k, double xmax); // uses gamma_dev1 or gamma_dev2, depending on k relative to LAMBDAMAX
double beta_variate(int n, int m, double xmin, double xmax);
double beta_variate1(int n, int m, double xmin, double xmax);
double beta_variate_0half(int n, int m);
double sample_1_minus_x_to_n(int n);
