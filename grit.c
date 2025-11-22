//  grit.c   main program
#include <mcheck.h>
#include <getopt.h>

#include "grit.h"
#include "params.h"
#include "structs.h"
#include "path_etc.h"
#include "path_tree_etc.h"
#include "histogram.h"
#include "perm.h"
#include "perm_etc.h" 
#include "chain.h"
#include "run_chains.h"
#include "signal.h"
#include "rngs.h"
#include "extdefs.h"

Permutation* simulate_rearrangement(Permutation* p1, long N_sim_revs, double sim_lambda_ratio);

int
main(int argc, char *argv[])
// main()
{
  long t_mainstart = time(NULL); // seconds from the epoch
  setvbuf(stdout, NULL, _IOLBF, BUFSIZ); // not sure what this does






  //*******************************************************


 // ***** process command line *****
  /* if (argc < 2) { */
  /*   fprintf(stderr, "Usage:  %s -in <dosages_file> [-out <output_filename> -nhgmr_max <max_nhgmr>] \n", argv[0]); */
  /*   print_usage_info(stderr); */
  /*   fprintf(stderr, "%d\n", (int)EXIT_FAILURE); */
  /*   exit(EXIT_FAILURE); */
  /* } */

  char* cl_genomes_filename = NULL;
  //FILE *g_stream = NULL; // for reading in genotypes
  char* cl_control_filename = NULL;
  //FILE *c_stream = NULL; // for reading in pedigrees (if specified)
  char* output_prefix = "";
  long cl_rng_seed = 0; // if left at this value, use value in control file.
  // static int cl_signed_flag = -1; // 1: signed, 0: unsigned, other: use control file value
  // static int cl_unsigned_flag = 0; // if set, unsigned, if neither set, use value in control file
  int cl_choose_signs = -1;
  
  int c;
  while(1){
    int option_index = 0;
    static struct option long_options[] = {
      {"control", required_argument, 0, 'c'}, // filename of control file 
      {"data",   required_argument, 0,  'd'}, // filename of genomes data set   
      {"output_prefix",  required_argument, 0,  'o'}, // output filename
      {"prefix", required_argument, 0, 'o'}, 
      {"seed", required_argument, 0, 'r'},
      {"rng_seed", required_argument, 0, 'r'},
      {"signed", no_argument, 0, 's'}, // --signed  set Choose_signs = 0  (i.e. signed)
      {"unsigned", no_argument, 0, 'u'}, // --unsigned
      {0,         0,                 0,  0 } // so will abort if unrecognized option
    };
     
    // c = getopt_long_only(argc, argv, "", long_options, &option_index);
    c = getopt_long(argc, argv, "", long_options, &option_index);
    printf("c: %i\n", c);
    if(c == -1) break;
    switch(c){
    case 'c':
      cl_control_filename = optarg;
      break;
    case 'd':
      cl_genomes_filename = optarg;
      break;
    case 'o':
      output_prefix = optarg;
    /* case 'r': */
    /*   cl_rng_seed = atoi(optarg); */
    /*   break; */
    case 'r':
      cl_rng_seed = atoi(optarg);
      break;
    case 's':
      cl_choose_signs = 0; // 
      break;
    case 'u':
      cl_choose_signs = 3;
      break;
    case '?':
      fprintf(stderr, "? case in command line processing switch. Exiting\n");
      exit(EXIT_FAILURE);
    default:
      fprintf(stderr, "default case. Exiting.\n");
      exit(EXIT_FAILURE);
    }
  }

  printf("cl control file: %s\n", cl_control_filename);
  printf("cl genomes file: %s\n", cl_genomes_filename);
  printf("output filename prefix: %s\n", output_prefix);
  printf("cl rng seed %ld\n", cl_rng_seed);
  printf("cl choose signs: %i\n", cl_choose_signs);
  //   printf("cl unsigned flag: %i\n", cl_unsigned_flag);
  //getchar();

//  *********************************************************

  char* control_filename = (cl_control_filename == NULL)? "grit_control" : cl_control_filename;
    
//  ***************  read in parameters from control file  *****************
  FILE* fp_control = fopen(control_filename, "r");
  Run_info_in r_in;
  if (fp_control != NULL){
    printf("before input_run_parms\n");
    input_run_params(fp_control, &r_in);
    //  ********   override control file values with command line values  *******
    if(cl_genomes_filename != NULL) r_in.Genome_data_file = cl_genomes_filename;
    if(cl_rng_seed != 0) r_in.Rand_seed = cl_rng_seed;
    if(cl_choose_signs >= 0){ r_in.Choose_signs = cl_choose_signs; printf("r_in.Choose_signs: %d\n", (int)r_in.Choose_signs); }
    if(r_in.Rand_seed < 0){
      r_in.Rand_seed = time(NULL);  // time returns secs since epoch
      printf("Control file has negative value for Rand_seed.\n");
      printf("Using rng seed from clock:  %ld\n", r_in.Rand_seed);
    }
    if(r_in.MC3_max_pathlength > MAXPATHLENGTH){ r_in.MC3_max_pathlength = MAXPATHLENGTH; }
    check_run_params(&r_in);
    printf("seed: %ld\n", r_in.Rand_seed); //getchar();
    seedrand(r_in.Rand_seed); 
  
  } else{ // no control file - fail.
    printf("Couldn't open input file %s; exiting.\n", control_filename);
    exit(EXIT_FAILURE);
  }
  fclose(fp_control);

  output_run_params(stdout, &r_in); // getchar();
   
  printf("Genomes data file: %s\n", r_in.Genome_data_file);
    
//    ****************   Read in genome data from file  *******************
  printf("Genome input file: %s \n", r_in.Genome_data_file); 
  FILE* fp_data = fopen(r_in.Genome_data_file, "r");
  if(fp_data == NULL){ printf("Couldn't open data file %s for reading. exiting\n", r_in.Genome_data_file); exit(EXIT_FAILURE); }
  Permutation* genomes[2];
  int distsinfile = get_genomes_from_file(fp_data, genomes, FALSE);
 
  printf("after get_genomes_from_file \n"); 
  printf("distsinfile: %i \n", distsinfile);
  // r_in.Use_distances = (r_in.Use_distances && distsinfile);
  if(!r_in.Use_distances){ // set L_g fields of permutation to N+M, even if distances in file
    genomes[0]->L_g = genomes[0]->n_mark + genomes[0]->n_chrom;
    genomes[1]->L_g = genomes[1]->n_mark + genomes[1]->n_chrom;
  }

  Permutation* p1 = genomes[0];
  Permutation* p2 = (r_in.N_rev_fake > 0)? // simulate p2 from p1 if requested
    simulate_rearrangement(p1, r_in.N_rev_fake, r_in.lambda_ratio_fake) :
    genomes[1];

  print_perm(stdout, p1);
  print_perm(stdout, p2);
  print_oxford_grid(stdout, p1, p2); 
  printf("Use_distances: %i \n", r_in.Use_distances);
  printf("Genomes from file %s: \n",r_in.Genome_data_file);
  printf("Genome 1: \n"); print_perm_brief(stdout, p1); print_perm(stdout, p1);
  printf("Genome 2: \n"); print_perm_brief(stdout, p2); print_perm(stdout, p2);
  printf("lambda_ratio_fake: %g \n", r_in.lambda_ratio_fake);
  
  printf("p1, p2, ninv, ntrans: %i %i     %i %i \n", p1->n_inv, p1->n_trans, p2->n_inv, p2->n_trans); //getchar(); 
  print_oxford_grid(stdout, p1, p2);  // getchar();
  // get Lambda_max from lambda_max
  long NpM = p1->n_mark + p1->n_chrom;
  if(r_in.Update_lambdas_mode <= 1){
    r_in.Lambda_max = (double)(NpM*(NpM-1)/2)*r_in.lambda_max;
  }
  else{
    r_in.Lambda_max = get_Lambda(NpM, r_in.lambda_max, r_in.lambda_max);
  }
 
  // open output files
  char rawfilename[200];
  for(long i=0; i<r_in.N_temperatures; i++){
    // open files Raw1.out, Raw2.out etc. for Raw output from different T chains
    sprintf(rawfilename, "%sT%iraw.out", output_prefix, (int)i);
    fp_raw[i] = fopen(rawfilename, "w");
    output_run_params(fp_raw[i], &r_in);
  }
  fflush(0);
        
  printf("****************************\n");
  print_perm(stdout, p1);
  print_perm(stdout, p2); 
  Cycle_decomposition a_cd;
  get_CD(p1, p2, &a_cd); printf("N, M, c_i, c_g: %i %i %i %i \n", a_cd.n_mark, a_cd.n_chrom, a_cd.n_int_cycles, a_cd.c_g);
  get_CD(p2, p1, &a_cd); printf("N, M, c_i, c_g: %i %i %i %i \n", a_cd.n_mark, a_cd.n_chrom, a_cd.n_int_cycles, a_cd.c_g);
     
  Chain_set* the_chain_set = construct_chain_set(p1, p2, &r_in);	     

  printf("*******************************************************************\n");
  printf("**********In main. chain_set constructed; now run chains **********\n");
  printf("*******************************************************************\n");

  long t_conv;
  long n_conv = run_chains(the_chain_set, &r_in, &t_conv);

  long t_mainstop = time(NULL);
  fprintf(stdout, "CPU time to convergence: %ld seconds\n", (long)t_conv);
  printf("MC updates to convergence:  %ld \n", n_conv);
  printf("number of path updates with log_p_a not(normal or zero): %i \n", n_path_p_a_nnoz);
  printf("number of lambda updates with log_p_a not(normal or zero): %i \n", n_lambdaIT_p_a_nnoz);
  printf("number of r/xi updates with log_p_a not(normal or zero): %i \n", n_rxi_p_a_nnoz);
  printf("Total CPU time (seconds): %ld \n\n", t_mainstop - t_mainstart);
  printf("Clocks, calls, clocks/call: \n");
  printf("Run_around: %g %g %g \n", run_around_clocks, run_around_calls, run_around_clocks/run_around_calls);
  printf("Run_around_old: %g %g %g \n", run_around_old_clocks, run_around_old_calls, run_around_old_clocks/run_around_old_calls);
  fflush(NULL);
  permfree(&p1); permfree(&p2);
  
  check_alloc_info(stdout);

  print_other_output(); //

  WAITFORKEY; // does nothing in linux case
  return 0;
}  // end of main

void check_alloc_info(FILE* fp)
{
  if(n_perm_alloc == 0 && n_p_alloc == 0 && n_step_alloc == 0 && n_path_alloc == 0 &&
     n_path_tree_alloc == 0 && n_path_tree_node_alloc == 0 && n_revseq_alloc == 0){
    fprintf(fp, "Checked that allocated memory freed, all ok. \n");
  }
  else{        
    fprintf(fp, "n perm alloc: %li \n",
	    n_perm_alloc);
    fprintf(fp, "n p alloc: %li \n",
	    n_p_alloc);
    fprintf(fp, "n step alloc: %li \n",
	    n_step_alloc);
    fprintf(fp, "n path alloc: %li \n",
	    n_path_alloc);
    fprintf(fp, "n path tree alloc: %li \n",
	    n_path_tree_alloc);
    fprintf(fp, "n path tree node alloc: %li \n",
	    n_path_tree_node_alloc);
    fprintf(fp, "n revseq alloc: %li \n",
	    n_revseq_alloc);
  }
} // end of function check_alloc_info




void handle_control_c(int sig)
{
  // set flag sigint_raised, so can stop soon (at convenient spot)
  printf("SIGINT raised; will stop soonish\n");
  sigint_raised = 1;
}



int
compare_ints(const void* aa, const void* bb)
{
  int a = *(int *)aa;
  int b = *(int *)bb;
  if(a > b) return 1;
  else if(a < b) return -1;
  else return 0;
}

int
compare_doubles(const void* aa, const void* bb)
{
  double a = *(double *)aa;
  double b = *(double *)bb;
  if(a > b) return 1;
  else if(a < b) return -1;
  else return 0;
}


double
trimean(const int* array, int n)
{
  // assumes an array in order
  double q1, q2, q3;

  q1 = (double)(array[n/4] + array[(n-1)/4])/2.0;
  q3 = (double)(array[(3*n/4)] + array[(3*n-1)/4])/2.0;
  q2 = (double)(array[(n-1)/2] + array[n/2])/2.0;
    
  return 0.25*(q1 + 2.0*q2 + q3);
}

double
rmean(int* Nconv, int n, double* sigma, double f)
{
  // throw out points > f*sigma from mean
  // recalculate mean, variance
  double sum = 0.0, sumsq = 0.0;
  double mean, result, v, rv;
  int i, nr=0;
  printf("top of rmean\n");
  for(i=0; i<n; i++){
    sum += (double) Nconv[i];
    sumsq += (double)Nconv[i]*(double)Nconv[i];
  }

  mean = sum/(double)n;
  v = (sumsq - mean*sum)/(double)(n-1);
  printf("middle of rmean, mean, v: %g %g \n",  mean, v);
  sum = 0.0; sumsq = 0.0; 
  for(i=0; i<n; i++){
    if(fabs((double)Nconv[i] - mean) < f*sqrt(v))
      {
	nr++;
	sum += (double) Nconv[i];
	sumsq += (double)Nconv[i]*(double)Nconv[i];
      }
  }
  result = sum/(double)nr;
  rv = (sumsq - result*sum)/(double)(nr-1);
  *sigma = sqrt(rv/(double)nr);
  printf("bottom of rmean, nr, rmean, rstddevmean: %i %g %g \n", nr, result, *sigma);
  return result;
}  

void print_other_output(void){
  // print out some other stuff - put into separate files, separate subroutines
  FILE* fpx;
  int cume=0;
  fpx = fopen("acc_etc.out", "w");

  fprintf(fpx, "# zeta_hist (dist of where trans are along path): \n");
  for(int ijk=0; ijk<20; ijk++){
    fprintf(fpx, "%g %i %i \n", 0.05*(ijk+0.5), zeta_hist_a[ijk], zeta_hist_b[ijk]);
    //    fprintf(stdout, "%g %i %i \n", 0.05*(ijk+0.5), zeta_hist_a[ijk], zeta_hist_b[ijk]);
  }fprintf(fpx, "\n");
        
  fprintf(fpx, "# nstucks: \n");
  for(int i=0; i<50; i++){
    fprintf(fpx, "%i  %i \n", i, nstucks[i]);
  }fprintf(fpx, "\n");
        
  fprintf(fpx, "# prop_lengths: \n");        
  for(int ijk = 0; ijk<2*MAXPATHLENGTH; ijk++){
    cume += prop_lengths[ijk];
    fprintf(fpx, "%5i %5i %5i \n", ijk, prop_lengths[ijk], cume);
  }
  fprintf(fpx, "\n");

  fprintf(fpx, "# lpa_hist: \n");        
  for(int ijk = 0; ijk<100; ijk++){
    fprintf(fpx, "%5i %g %i \n", ijk, 0.5*(ijk - 70 + 0.5), lpa_hist[ijk]);
  } fprintf(fpx, "\n");

        
   

  fprintf(fpx, "#track length i, ptl[i], atl[i]\n");
  {
    double sum_ptl = 0.0, sum_atl = 0.0;
    for(int ijk=0; ijk<MAX_N_MARKERS_PER_CHROMOSOME; ijk++){
      sum_ptl += ptl[ijk]; sum_atl += atl[ijk];
    }
    for(int ijk=0; ijk<40; ijk++){
      fprintf(fpx, "%i %g %g \n", ijk, ptl[ijk]/sum_ptl, atl[ijk]/sum_atl);
    }
  }
}

// generate p2 from p1 by random reversals
Permutation* simulate_rearrangement(Permutation* p1, long N_sim_revs, double sim_lambda_ratio){ 
  Reversal the_rev;
  Cycle_decomposition the_cd;
  int li = 0, lt = 0;
  double prob;    	
  Permutation* p2 = copy_perm(p1);
  for (long i=0; i<N_sim_revs; i++){
    get_CD(p2, p2, &the_cd);
    the_rev = rand_reverse_r(p2, sim_lambda_ratio, &prob, the_cd.elements);
    if(the_rev.is_inversion == TRUE) li++;
    else lt++;
    // printf("the_rev.left, right: %g %g \n", the_rev.ld, the_rev.rd);
  }
  //print_perm_brief(stdout, p1); print_perm_brief(stdout, p2);
  printf("Inversions, translocations in simulated data:  li, lt: %i %i \n", li, lt);   // getchar();
  return p2;
}
