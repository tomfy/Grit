// histogram.c    Histogram stuff

#include "grit.h"
#include "structs.h"
#include "exts.h"
#include "path.h"
#include "histogram.h"
#include "gr_misc.h"


Histogram* 
create_histogram(int nbins, double minval, double binwidth)
{
  Histogram* hist;

  hist = (Histogram*)chmalloc(sizeof(Histogram));
  hist->nbins = nbins;
  hist->minval = minval; // lower edge of leftmost bin, i.e. smallest value that goes in bin hist->histogram[0]
  hist->binwidth = binwidth;
  hist->sum = 0.0;
  hist->peak = 0.0;
  hist->peak_bin = -1;
  hist->mode = UNKNOWN;
  hist->low_ne = nbins+10;
  hist->hi_ne = -1; 
  hist->histogram = (double*)chcalloc((size_t)nbins, sizeof(double));
  hist->stats = (Stats*)chmalloc(sizeof(Stats));
  initialize_stats(hist->stats);
  return hist;
}

void
cred_int(Histogram* hist, double CL, double* xlo, double* xhi, double* ylo, double* yhi)
{
    int i=0,j=hist->nbins-1;
    double* a=hist->histogram; 
    double prob=0.0, newprob, pi, pj;
    double sum = hist->sum, binwidth = hist->binwidth, tailprob=1.0-CL;

    *xlo = hist->minval; *xhi = *xlo + binwidth*hist->nbins;
    *ylo = *yhi = 0.0;
    
    do{
        pi = a[i]/sum;
        pj = a[j]/sum;
        if(pi < 0.0 || pj < 0.0) {
            printf("in cred_int, negative probability\n");
        }
        if (pi < pj){
            newprob = prob + pi;
            *ylo = pi/binwidth;
            if (newprob > tailprob){ // *xlo will lie in ith bin
                *xlo = *xlo + (tailprob - prob)/pi*binwidth; // printf("xlo, xhi: %g %g \n", *xlo, *xhi);
                    // *ylo = pi/binwidth; *yhi = pj/binwidth;
                return;
            }
            
            prob = newprob;
            i++; *xlo += binwidth;
        }
        else if (pj <= pi){
            newprob = prob + pj;
            *yhi = pj/binwidth;
            if (newprob > tailprob){ // *xhi will lie in jth bin
                *xhi = *xhi - (tailprob - prob)/pj*binwidth; // printf("xlo, xhi: %g %g \n", *xlo, *xhi);
                    // *ylo = pi/binwidth; *yhi = pj/binwidth;
                return;
            }
            prob = newprob;
            j--; *xhi -= binwidth; // printf("xhi: %g \n", *xhi);
        }
    }while(TRUE);
}  // end of cred_int

void
cred_int1(double* x, double* a, int N, double CL, double* xlo, double* xhi, double* ylo, double* yhi)
{
        // gets credible interval with confidence level CL, given arrays x and a
    int i,j;
    double prob=0.0, newprob, pi, pj;
    double sum = 0.0, binwidth = x[1] - x[0], tailprob=1.0-CL;

    *xlo = x[0] - 0.5*binwidth; *xhi = x[N-1] + 0.5*binwidth;
    *ylo = *yhi = 0.0;

    for(i=0; i<N; i++){
        sum += a[i];
    }
        
    i=0; j=N-1;
    do{
        pi = a[i]/sum;
        pj = a[j]/sum;
        if(pi < 0.0 || pj < 0.0) {
            printf("in cred_int, negative probability\n");
        }
        if (pi < pj){
            newprob = prob + pi;
            *ylo = pi/binwidth;
            if (newprob > tailprob){ // *xlo will lie in ith bin
                *xlo = *xlo + (tailprob - prob)/pi*binwidth; printf("xlo, xhi: %g %g \n", *xlo, *xhi);
                    // *ylo = pi/binwidth; *yhi = pj/binwidth;
                return;
            }
            
            prob = newprob;
            i++; *xlo += binwidth;
        }
        else if (pj <= pi){
            newprob = prob + pj;
            *yhi = pj/binwidth;
            if (newprob > tailprob){ // *xhi will lie in jth bin
                *xhi = *xhi - (tailprob - prob)/pj*binwidth; printf("xlo, xhi: %g %g \n", *xlo, *xhi);
                    // *ylo = pi/binwidth; *yhi = pj/binwidth;
                return;
            }
            prob = newprob;
            j--; *xhi -= binwidth; // printf("xhi: %g \n", *xhi);
        }
    }while(TRUE);
}  // end of cred_int1
        
    

void
free_histogram(Histogram* hist)
{
    free(hist->stats);
  free(hist->histogram);
  free(hist);
}

void 
insert_in_histogram(Histogram* hist, double value, double weight)
{
  int bin = (int)((value - hist->minval)/(hist->binwidth));
  if (bin < 0) bin = 0;
  if (bin >= hist->nbins) bin = (hist->nbins)-1;
  hist->histogram[bin] += weight;
  hist->sum += weight;
  if(weight > 0.0){
      if(bin > hist->hi_ne) hist->hi_ne = bin;
      if(bin < hist->low_ne) hist->low_ne = bin;
  }
  if(hist->histogram[bin] > hist->peak){
      hist->peak = hist->histogram[bin];
      hist->peak_bin = bin;
      hist->mode = hist->minval + (0.5 + bin)*hist->binwidth;
  }
//  printf("in insert_in_hist. before insert_in_stats. value, weight: %g %g \n", value, weight);
  insert_in_stats(hist->stats, value, weight);
}

void
print_histogram(FILE* fp, Histogram* hist, int normalize, double A)
{ // print a histogram
        // normalize = 0    don't normalize
        // normalize = 1  normalize so sum of bins is A
        // normalize = 2   normalize so peak is A
        // normalize = 3  normalize so sum of bins is A/binwidth (i.e. integral is A)
    int i;
    double denom = 1.0;
    double v, max = UNKNOWN;
    double* h;
    int lo_print_bin = 0, hi_print_bin = hist->nbins - 2;


    if(hist->sum > 0.0){
        for(i=1; i<hist->nbins-1; i++){
            if(hist->histogram[i] > 0){ lo_print_bin = i-1; break; }
        }
        for(i=hist->nbins-3; i>0; i--){
            if(hist->histogram[i] > 0){ hi_print_bin = i+1; break; }
        }
        
        if(normalize == 1){
            if(hist->sum > 0.0)denom = hist->sum/A;
            else {
                printf("In print_histogram. Can't normalize histogram, sum is <= 0\n");
                denom = 1.0;
            }
        }
        else if(normalize == 2){
            for(i=0, h = hist->histogram; i<hist->nbins; i++, h++){
                v = *h;
                if (v > max) max = v;
            }
            if(max > 0.0){
                denom = max/A;
            }
            else  {
                printf("In print_histogram. Can't normalize histogram, max is <= 0\n");
                denom = 1.0;
            }
        }
        else if(normalize == 3){
            if(hist->sum > 0.0)denom = hist->sum*hist->binwidth/A;
            else {
                printf("In print_histogram. Can't normalize histogram, sum is <= 0\n");
                denom = 1.0;
            }
        }
      
    }
    else{
        lo_print_bin = 0; hi_print_bin = 1; denom = 1.0;
    }

  
 
  
    for(i=lo_print_bin; i<=hi_print_bin; i++){ // last bin is overflow - don't print
        fprintf(fp, "%g %g %g\n", hist->minval + hist->binwidth*(i + 0.5), hist->histogram[i], hist->histogram[i]/denom);
    }
    fprintf(fp, "\n");
}

void
print_histograms(FILE* fp, int Nhists, Histogram** the_hists, int normalize, double A)
{ // print a histogram
        // normalize = 0    don't normalize
        // normalize = 1  normalize so sum of bins is A
        // normalize = 2   normalize so peak is A
        // normalize = 3 normalize so integral is A
    int i, j;
    double v, max = UNKNOWN;
    double* h;
    Stats the_stats;
    Histogram* hist;
    double denom[NCHAINMAX];
    int lo_print_bin[NCHAINMAX], hi_print_bin[NCHAINMAX];
    int lpb, hpb;
    int non_empty = FALSE; // 

    for(j=0; j<Nhists; j++){
        if(the_hists[j]->sum > 0.0){ non_empty = TRUE; }
    }
  
    if(non_empty){
        for(j=0; j<Nhists; j++){ 
            for(i=1; i<the_hists[j]->nbins-1; i++){
                if(the_hists[j]->histogram[i] > 0){ lo_print_bin[j] = i-1; break; }
            }
            for(i=the_hists[j]->nbins-3; i>0; i--){
                if(the_hists[j]->histogram[i] > 0){ hi_print_bin[j] = i+1; break; }
            }
        }
        lpb = lo_print_bin[0]; hpb = hi_print_bin[0];
        for(j=1; j<Nhists; j++){
            if(lo_print_bin[j] < lpb) lpb = lo_print_bin[j];
            if(hi_print_bin[j] > hpb) hpb = hi_print_bin[j];
        }

      
            // normalize all the histograms
        for (j=0; j<Nhists; j++){
            hist = the_hists[j];
            if(normalize == 1){
                if(hist->sum > 0.0)denom[j] = hist->sum/A;
                else {
                    printf("In print_histogram. Can't normalize histogram, sum is <= 0\n");
                    denom[j] = 1.0;
                }
            }
            else if(normalize == 2){
                for(i=0, h = hist->histogram; i<hist->nbins; i++, h++){
                    v = *h;
                    if (v > max) max = v;
                }
                if(max > 0.0){
                    denom[j] = max/A;
                }
                else  {
                    printf("In print_histogram. Can't normalize histogram, max is <= 0\n");
                    denom[j] = 1.0;
                }
            }
            else if(normalize == 3){
                if(hist->sum > 0.0)denom[j] = hist->sum*hist->binwidth/A;
                else {
                    printf("In print_histogram. Can't normalize histogram, sum is <= 0\n");
                    denom[j] = 1.0;
                }
            }
        }
    }
    else{
        printf("In print_histograms. All histograms are empty. \n");
        lpb = 0; hpb = 1;
        for(j=0; j<Nhists; j++){ denom[j] = 1.0; }
    }

       
    for(i=lpb; i<=hpb; i++){
        initialize_stats(&the_stats);
        for(j=0; j<Nhists; j++){          
            insert_in_stats(&the_stats, the_hists[j]->histogram[i]/denom[j], 1.0);
        }
        fprintf(fp, "%g %g %g ", the_hists[0]->minval +  the_hists[0]->binwidth*(i + 0.5),
                the_stats.mean, the_stats.stddev_mean);
     
        for(j=0; j<Nhists; j++){
            fprintf(fp, "%g ", the_hists[j]->histogram[i]/denom[j]);
        }fprintf(fp, "\n");
    }fprintf(fp, "\n");
//  fprintf(fp, "%g %g \n\n", the_hists[0]->minval + the_hists[0]->binwidth*(20.0*hist->nbins + 0.5), 0.0);
} // end print_histograms

void
initialize_stats(Stats* stat)
{
    stat->N = 0; stat->sum = 0.0; stat->sumsqr = 0.0;
}

void
insert_in_stats(Stats* stat, double value, double weight)
{
    double mean;
    stat->N += weight;
    stat->sum += weight*value;
    stat->sumsqr += weight*value*value;
    mean = stat->sum/stat->N;
    stat->mean = mean;
        //  printf("mean: %g \n", mean); getchar();
    if(stat->N>1){
    stat->stddev_mean = sqrt((stat->sumsqr/stat->N - mean*mean)/(stat->N-1));
    stat->stddev = stat->stddev_mean*sqrt(stat->N);
    }
    else{
        stat->stddev = stat->stddev_mean = UNKNOWN;
    }
}

void print_stats(FILE* fp, Stats* stats)
{
    fprintf(fp, "N: %g; mean: %g; sigma: %g; sigma_mean: %g \n", stats->N, stats->mean, stats->stddev, stats->stddev_mean);
}

