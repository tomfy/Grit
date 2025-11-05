// functions related to gamma function, etc.
// use lambda = gamma_dev(the_path->length + 1, Lambda_max) to update lambda

//loggamma function, etc

#include "grit.h"
#include "gamma_dev.h"
#include "structs.h"
#include "exts.h"
#include "rngs.h"

double loggamma(double xx) // based on Numerical Recipes routine
{
    double cof[6] = {76.18009173, -86.50532033, 24.01409822,
                     -1.231739516, 0.120858003e-2, -0.536382e-5};
    double stp = 2.50662827465;
    double x, tmp, ser;
    int i;

    x = xx - 1.0;
    tmp = 5.5+x;
    tmp = (x + 0.5)*log(tmp) - tmp;
    ser = 1.0;
    for(i=0; i<6; i++){
        x += 1.0;
        ser += cof[i]/x;
    }
    tmp += log(stp*ser);
    
    return tmp;
}

double gseries(double a, double x) // based on Numerical Recipes routine
{
    int itmax = 1000, n;
    double gln = loggamma(a);
    double ap, sum, del, eps = 1.0e-12;
    
    //printf("in gseries; a, x: %f %f \n", a, x);
    
    if(x < 0.0){
        printf("in gseries, x < 0.0 \n"); exit(EXIT_FAILURE);
    }
    else if (x == 0.0){
        return 0.0;
    }
    else {
        ap = a;
        sum = 1.0/a;
        del = sum;
        for (n=0; n<itmax; n++){
            ap += 1.0;
            del *= (x/ap);
            sum += del;
            if(fabs(del) < (fabs(sum)*eps)) return sum*exp(-x + a*log(x) - gln);
        }
        printf("in gseries; No convergence after: %i iterations. a: %g, x: %g, test: %g \n", itmax, a, x, fabs(del/sum));
        exit(EXIT_FAILURE);
    }
}

double gcf(double a, double x) // based on Numerical Recipes routine
{ 
    int itmax = 1000;
    double eps = 1.0e-12;
    double gln = loggamma(a);
    double gold = 0.0;
    double a0 = 1.0;
    double a1 = x;
    double b0 = 0.0;
    double b1 = 1.0;
    double fac = 1.0;
    double n, ana, anf, g;
    double test;
    static int first = TRUE;

    if(first){
        printf("gcf called for first time: itmax, a, x: %i %g %g\n", itmax, a, x);
        first = FALSE;
    }
    for (n=1.0; (int)n <= itmax; n += 1.0){
        ana = n - a;
        a0 = (a1 + a0*ana)*fac;
        b0 = (b1 + b0*ana)*fac;
        anf = n*fac;
           
        a1 = x*a0 + anf*a1; 
        b1 = x*b0 + anf*b1;
          
        if(a1 != 0.0){
            fac = 1.0/a1;
            g = b1*fac;
               
            test = (g - gold)/g;
            if(fabs(test) < eps) {
                return exp(-x + a*log(x) - gln)*g;
            }
            gold = g;
        }
        else printf("In gcf, a1 = %g \n", a1);
        
    }    
      printf("in gcf; No convergence after %i iterations; a: %g x: %g test: %g\n", itmax, a, x, test);
	exit(EXIT_FAILURE);
}

double incomplete_gamma(double a, double x)    // based on Numerical Recipes routine
{
        // incomplete_gamma(a,x) = 1/(a-1)! * integral(e^-t * t^(a-1) dt) from 0 to x
    
    if((x < 0.0) || (a <= 0.0)) printf("in incomplete_gamma, x<0 or a<=0\n");
    
    if (x < a+1.0) {
        return gseries(a, x);
    }
    else{
        //printf("In incomplete gamma, gcf branch. a, x: %g %g \n", a, x);
        return 1.0 - gcf(a, x);
    }
}

double gamma_dev2(double a, double xmax) // value returned is < xmax
{ // uses Newton's method to find x such that incomplete_gamma(a, x) = y_targ
  // y_targ is random number between 0 and incomplete_gamma(a, xmax)
  int i = 0, itmax = 100;
  double eps = 1.0e-12;
  double x, y;
      //double y_targ = (double)RANDF()/((double)RANDMAX + 1)*incomplete_gamma(a, xmax);
  double y_targ = drand(&rng_long)*incomplete_gamma(a, xmax);

  x = xmax; // initial guess to feed into Newton's method0
  while(fabs(y = incomplete_gamma(a, x) - y_targ) > eps){
          //printf("in gamma_dev2, in loop. i: %i \n", i);
      
    x = x - y*exp(x - (a - 1.0)*log(x) + loggamma(a));
    i++; 
    if(i > itmax) {
      printf("in gamma_dev2, i > itmax (newtons' method) \n");
      exit(EXIT_FAILURE);
    }
  } 
  return x;
}

double gamma_dev1(int k, double xmax)

// based on routine in Numerical Recipes 
// for small k use sum of k exponentially distributed values, 
// for large k use rejection method with lorentzian as 
// comparison function 
// k must be positive

// returns a value x which samples the x < xmax part of the pdf 
// P(x) = exp(-x)*x^(k-1)/(k-1)! 

{ 
    int i, done; 
    double x = 0.0, y, km, s, e; 
 
    if (k < 1) { 
        printf("gamma_dev1 called with non-positive argument. \n"); 
        exit(EXIT_FAILURE); 
    } 
 
    do{ 
        if (k <= 7){
            x = 0.0; // need this in case x>xmax and loop repeats!
            for (i=0; i< k; i++){
                x -= log(drand(&rng_long)); 
            } 
        } 
        else{ 
            do{ 
                done = FALSE;
                   
                x = drand(&rng_long); 
                if (x > 0.0){ 
                    y = tan(M_PI*(x-0.5)); 
                    km = (double)(k-1); 
                    s = sqrt(2.0*km + 1); 
                    x = s*y + km; 
                    if (x > 0.0){ 
                        e = (1.0 + y*y)*exp(km*log(x/km) - s*y);
                        done = drand(&rng_long) <= e; 
                    } 
                }
            } 
            while (!done); 
        } 
    } 
    while (x > xmax); 
    return x; 
} 

double gamma_dev(int k, double xmax)
{
  double a = (double)k;

  if(xmax < a - 1.4*sqrt(a)){
    return gamma_dev2(a, xmax);
  }
  else{
    return gamma_dev1(k, xmax);
  }
}

double sample_1_minus_x_to_n(int n)
{
        // returns a sample from p(x) = (1-x)^n    (0 <= x <= 1/2)
    double x;

    while(1==1){
        x = pow(drand(&rng_long), 1.0/((double)n + 1.0));
        if(x >= 0.5) return 1.0 - x;
    }
}


double beta_variate_0half(int n, int m)
{
        //  samples beta distribution, in range 0 to 1/2
        // x^n*(1-x)^m
    
        //   printf("top of beta_variate_0half. n, m: %i %i \n", n, m);
    if(n>m){
        if(m==0) return 0.5*pow(drand(&rng_long), 1.0/(double)(n+1)); // highly unusual branch
        else return beta_variate1(n, m, 0.0, 0.5); }
    else if(n == 0){
        return sample_1_minus_x_to_n(m);} 
    else{
        return beta_variate(n, m, 0.0, 0.5);
    }
        //  printf("bottom of beta_variate_0half\n");
}

double
beta_variate(int n, int m, double xmin, double xmax)
{
        // returns a sample from probability distribution with
        // (unnormalized) density p(x) = x^n*(1-x)^m  for xmin <= x <= xmax
        // 0 <= xmin <= xmax <= 1.0
        // requires many iterations if n/(n+m) (i.e. xhat) is much outside range defined by xmin, xmax
        // doesn't work for n = 0 or for m = 0.
    int N = n+m;
    double rn = (double)n;
    double rm = (double)m;
    double xhat = rn/(double)N;
    double p0 = pow(xhat, rn)*pow((1.0-xhat), rm);
    double logp0 = rn*log(xhat) + rm*log(1.0 - xhat);
    double alpha = (double)N*sqrt(1.0 - 0.5*((rn-1.0)/rn + (rm-1.0)/rm));
    double a0 = atan(alpha*(xmin - xhat)), a1 = atan(alpha*(xmax - xhat));
        //, d = a1 - a0;
    double yr;
    double ff = 0.98; //
    double x;
    double q, p;
    double random;
    int itcount = 0;
    int itmax = 1000;

        // printf("top of beta_variate1. n, m, xmin, xmax: %i %i %g %g \n", n, m, xmin, xmax);

          if((double)n/(double)(n+m) > xmax) {printf("in beta_variate, xhat > xmax; use beta_variate1 instead\n"); getchar();}
        // x = tan(drand(&rng_long)*a1 + (1.0-yr)*a0)/alpha + xhat;

        //   if(xmin < 0.0 || xmax < xmin || xmax > 1.0)
    while(TRUE){
        itcount++;
        yr = drand(&rng_long);
        x = tan(yr*a1 + (1.0-yr)*a0)/alpha + xhat;
        q = 1.0/(1.0 + pow(alpha*(xhat - x), 2.0));
        p = ff*pow(x, rn)*pow(1.0 - x, rm)/p0;
        double logp = log(ff) + rn*log(x) + rm*log(1.0 - x) - logp0;
      
        random = drand(&rng_long);
            //   printf("%g %g %g %g %g  %i \n", q, p, p/q, q-p, random, (random < p/q)); 
            //
            //  if(random < p/q){ // accept
        if(log(q*random) < logp){ // accept
                //  printf("itcount: %i \n", itcount);
                //   printf("%i\n", itcount);
            return x;
        }
        if (itcount > itmax){
            printf("in beta_variate. itcount too big. itcount = %i ,n, m: %i %i \n", itcount, n, m);
                // cerr << "in beta_variate, itcount: " << itcount << "  too big" << endl;
            getchar();
        }
    }
        //  printf("bottom of beta_variate1\n");

}

double
beta_variate1(int n, int m, double xmin, double xmax)
{
        // returns a sample from probability distribution with
        // (unnormalized) density p(x) = x^n*(1-x)^m  for xmin <= x <= xmax

        // use stepwise constant function as comparison function
        // 0 <= xmin <= xmax <= 1.0

    double sigma = sqrt((double)(n*m))/pow((double)(n+m), 1.5);
    double step = 0.5*sigma;

    double x[100];
    double f[100];
    double p[100];

    double the_x = xmax, next_x;
    double sump = 0; 
    int i, j, nintervals, maxits = 100;
    double randn;
    int done = FALSE;
    
    if((double)n/(double)(n+m) < xmax) {
        printf("In beta_variate1, xhat is < xmax\n"); getchar();
    }
   
    for(i=0, p[0] = 0.0; !done; i++){
             
        x[i] = the_x;
        f[i] = pow(the_x, (double)n)*pow(1.0 - the_x, (double)m);
        next_x = the_x - step;
        if(next_x <= xmin){
            next_x = xmin;
            done = TRUE;
        }          
        p[i+1] = p[i] + (the_x - next_x)*f[i];
            //   printf("i, p[i+1]:  %i %g \n", i, p[i+1]);
        
        the_x = next_x;
        if(i > 3) step *= 2.0;
    }
    nintervals = i;
    sump = p[i];
        //  getchar();

    for(j=0; j<maxits; j++){
        randn = drand(&rng_long)*sump;
    
        for(i=0; i<nintervals; i++){
            if(randn < p[i+1]){
                the_x = x[i] + (x[i+1]-x[i])*(randn - p[i])/(p[i+1] - p[i]);
                    //  printf("%g %g \n", randn, p[i+1]);
                if(drand(&rng_long)*f[i] < pow(the_x, (double)n)*pow(1.0 - the_x, (double)m)){
                        //   printf("beta_variate1 returning: %g \n", the_x);
                    return the_x;
                }
                else break;
            }
        }
    }
    printf("in beta_variate1, > %i iterations. Too many!\n", maxits);
}
            
