//  grit.c   main program
#include <mcheck.h>

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

// ########### streamline it by removing loops over param values and data sets, etc.
// 

#include "extdefs.h"

int main()
{
    int i, j, n_conv, NpM;
    double t_mainstart = 0.0, t_mainstop, cpu_time;
    Permutation* p1[NDATAMAX] = {NULL}, * p2[NDATAMAX] = {NULL}; // arrays of pointers to perms
    FILE* fp_in1, * fp_in2;
    Histogram* n_conv_histogram = create_histogram(1000, 0.0, 250.0);
    Histogram* t_conv_histogram = create_histogram(1000, 0.0, 0.25); 
   
    int Nconvdists[1][NDATAMAX] = {{0}}; // to store burn-in lengths for conv. studies
    Run_info_in r_in; // no longer global!
    double prob;
    Permutation* genomes[2];
    int Nconv[NDATAMAX] = {0};
    double tconv[NDATAMAX] = {0};
    int distsinfile;
    Chain_set* the_chain_set;

    t_mainstop = GETCPUTIME;
    
    setvbuf(stdout, NULL, _IOLBF, BUFSIZ);
    
/*     if(HANDLE_SIGINT)(void)signal(SIGINT, handle_control_c);  // signal returns previous handler (type signal_handler_t) in case we want to restore it.     */
/*     t_mainstart = GETCPUTIME; */

/*     if((fjunk = fopen("xxxx", "w")) != NULL){ */
/*         fprintf(fjunk, "ksflsfdlkfdkfdsklsdfklj\n"); */
/*     } */
/*     else */
/*         printf("couldn't open xxxx for writing\n"); */

/*     if((fpz = fopen("zzz", "r"))!= NULL){ */
/*         printf("opened zzz\n"); */
/*     } */
/*     else{ */
/*         printf("couldn't open zzz for reading\n"); */
/*      } */
        // read in parameters 
    fp_in1 = fopen("grit_control", "r");
    if (fp_in1 != NULL){
        printf("before input_run_parms\n");
        input_run_params(fp_in1, &r_in);
        printf("Using input parameters read from file.\n");
    }
    else{
        printf("Couldn't open input file grit_control; exiting.\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp_in1);

    
   /*      // read in genome data from file */
/*         // and if N_rev_fake > 0, generate fake data */

    printf("Genome input file: %s \n", r_in.Genome_data_file);
    fp_in2 = fopen(r_in.Genome_data_file, "r");
    distsinfile = get_genomes_from_file(fp_in2, genomes, FALSE);
        //  distsinfile = get_data_from_file_new(fp_in2, genomes, false);
        //     distsinfile = get_data_from_file_new(fp_in2, genomes, r_in.Use_distances);
    printf("after get_data_from_file_new... \n"); 
    printf("distsinfile: %i \n", distsinfile);
        // r_in.Use_distances = (r_in.Use_distances && distsinfile);
    if(!r_in.Use_distances){ // set L_g fields of permutation to N+M, even if distances in file
        genomes[0]->L_g = genomes[0]->n_mark + genomes[0]->n_chrom;
        genomes[1]->L_g = genomes[1]->n_mark + genomes[1]->n_chrom;
    }
    print_perm(stdout, genomes[0]);
    print_perm(stdout, genomes[1]);
    print_oxford_grid(stdout, genomes[0], genomes[1]);   // getchar();
    printf("Use_distances: %i \n", r_in.Use_distances); //getchar();
        //  getchar();
    printf("Genomes from file %s: \n",r_in.Genome_data_file);
    printf("Genome 1: \n"); print_perm(stdout, genomes[0]);
    printf("Genome 2: \n"); print_perm(stdout, genomes[1]);
    //getchar();
           //seed the random number generator
    if(r_in.Rand_seed < 0){
      r_in.Rand_seed = (unsigned long)time(NULL);
    }
    seedrand(r_in.Rand_seed);
    printf("rand seed: %li \n", r_in.Rand_seed);
    printf("lambda_ratio_fake: %g \n", r_in.lambda_ratio_fake);
    printf("drand(): %g \n", drand(&rng_long));
   
    output_run_params(stdout, &r_in);
    check_run_params(&r_in);
        //  getchar();
    
    {
        Reversal the_rev;
        Cycle_decomposition the_cd;
        int li = 0, lt = 0;
              int g1 = 0, g2 = 1;
            //   int g1 = 1, g2 = 0;
    for (j=0; j<r_in.N_data; j++){
        p1[j] = copy_perm(genomes[g1]);
        if(r_in.N_rev_fake >= 1){    
            p2[j] = copy_perm(genomes[g1]);
                // generate p2 from p1 by random reverses
            for (i=0; i<r_in.N_rev_fake; i++){
                get_CD(p2[j], p2[j], &the_cd);
                the_rev = rand_reverse_r(p2[j], r_in.lambda_ratio_fake, &prob, the_cd.elements);
                if(the_rev.is_inversion == TRUE) li++;
                else lt++;
                    //  printf("the_rev.left, right: %i %i \n", the_rev.left, the_rev.right);
                    //  print_perm(stdout, p2[j]);
            }
            print_perm(stdout, p1[j]); print_perm(stdout, p2[j]);
            printf("li, lt: %i %i \n", li, lt);
        }
        else p2[j] = copy_perm(genomes[g2]);
    }
    }
    permfree(genomes); permfree(genomes+1);
    printf("p1, p2, ninv, ntrans: %i %i     %i %i \n", p1[0]->n_inv, p1[0]->n_trans, p2[0]->n_inv, p2[0]->n_trans); //getchar(); 
        //  getchar();

    print_oxford_grid(stdout, p1[0], p2[0]); // getchar();
        // get Lambda_max from lambda_max
    NpM = p1[0]->n_mark + p1[0]->n_chrom;
    if(r_in.Update_lambdas_mode <= 1){
        r_in.Lambda_max = (double)(NpM*(NpM-1)/2)*r_in.lambda_max;
    }
    else{
        r_in.Lambda_max = get_Lambda(NpM, r_in.lambda_max, r_in.lambda_max);
    }

    {
        int ixx;
            //   printf("LTHEATKNEE, LTHEATFACTOR: %i %g \n", LTHEATKNEE, LTHEATFACTOR);
        for(ixx = 1; ixx<30; ixx++){
            printf("%i %g \n", ixx, f_Ltheat(r_in.T_hot, r_in.T_hot, ixx));
        }

            //   getchar();
    }

   

  
// open output files
    if(r_in.N_data == 1){

      /* fp_out5 = fopen("Chain_averages.out", "w"); // <L> for each chain */
      /* output_run_params(fp_out5, &r_in); */
        
      /* fp_out7 = fopen("Llambda_histos.out", "w"); */
      /* output_run_params(fp_out7, &r_in); */

      /* fp_LiLtout = fopen("LiLt.out", "w"); */
      /* output_run_params(fp_LiLtout, &r_in); */

      /* fp_lambdaIlambdaTout = fopen("lambdaIlambdaT.out", "w"); */
      /* output_run_params(fp_lambdaIlambdaTout, &r_in); */

        if(WRITERAW){
          /*   if(RUNCHOLD){ */
           
/*                 fp_out1 = fopen("Raw.out", "w"); */
/*                 output_run_params(fp_out1, &r_in); */
/*                 fp_rawhot = fopen("Raw_hot.out", "w"); */
/*                 output_run_params(fp_rawhot, &r_in); */
/*             } */
/*             else */
            {
                char rawfilename[200];
                for(i=0; i<r_in.N_temperatures; i++){
                        // open files Raw1.out, Raw2.out etc. for Raw output from different T chains
                    sprintf(rawfilename, "RawT%i.out", i);
                    fp_raw[i] = fopen(rawfilename, "w");
                    output_run_params(fp_raw[i], &r_in);
                }
            }
           // fp_rawq = fopen("Rawq.out", "w"); 
        }
    }
    
    /* fp_out4 = fopen("L_lambda_dist.out", "w"); // histograms of L and lambda. */
    /* output_run_params(fp_out4, &r_in);  */

    /* fp_out2 = fopen("XBW.out", "w"); // averages, W, B, etc. */
    /* output_run_params(fp_out2, &r_in); */
        
    /* fp_out3 = fopen("Converge_summary.out", "w"); // N_conv, t_conv etc. */
    /* output_run_params(fp_out3, &r_in); */

    /* fp_out8 = fopen("translocations.out", "w"); // N_conv, t_conv etc. */
    /* output_run_params(fp_out8, &r_in); */

    /* fp_out9 = fopen("Conv_Summary.out", "w"); */
    /* output_run_params(fp_out9, &r_in); */

    /* fp_out10= fopen("Conv_dists.out", "w"); */
    /* output_run_params(fp_out10, &r_in); */

    printf(" \n");
  
    fflush(0);

        // 
    if(r_in.MC3_max_pathlength > MAXPATHLENGTH){ r_in.MC3_max_pathlength = MAXPATHLENGTH; }

    seedrand(r_in.Rand_seed);  
        
        printf("****************************\n");
        print_perm(stdout, p1[0]);
        print_perm(stdout, p2[0]);
        {   
            Cycle_decomposition a_cd;
            get_CD(p1[0], p2[0], &a_cd);
            printf("N, M, c_i, c_g: %i %i %i %i \n", a_cd.n_mark, a_cd.n_chrom, a_cd.n_int_cycles, a_cd.c_g);
            get_CD(p2[0], p1[0], &a_cd);
            printf("N, M, c_i, c_g: %i %i %i %i \n", a_cd.n_mark, a_cd.n_chrom, a_cd.n_int_cycles, a_cd.c_g);
        }
       
            // run N_chains  chains to convergence
	j=0; // used to be able to analyze multiple data sets in a single run
	// and j indexed the data set.
     
            {               
                    //   printf("about to call run_chains_new. \n"); getchar();
	      print_perm(stdout, p1[j]);
	      print_perm(stdout, p2[j]);
              the_chain_set = construct_chain_set(p1[j], p2[j], &r_in);
                
    printf("**************************************************************\n");
    printf("*************In main. After construct_chain_set **************\n");
    printf("**************************************************************\n");
              
                n_conv = run_chains_new(the_chain_set, &r_in, &cpu_time);
            }
	    
             // "backward"
             //      n_conv = run_chains(p2[j], p1[j], &r_in, &cpu_time);
            Nconv[j] = n_conv;
            tconv[j] = (double)cpu_time;
            
            insert_in_histogram(n_conv_histogram, (double)n_conv, 1.0);
            insert_in_histogram(t_conv_histogram, (double)cpu_time, 1.0);
            fflush(NULL);
                //   printf("bottom of loop\n");
             
	    //  } // end of loop over N_data data sets
        qsort(Nconv, r_in.N_data, sizeof(int), compare_ints);
        qsort(tconv, r_in.N_data, sizeof(double), compare_doubles);
	for(j=0; j<r_in.N_data; j++){
	  printf("j, Nconv[j]: %i %i \n", j, Nconv[j]);
	  Nconvdists[0][j] = Nconv[j];
        }
	
         {
         double median, q1, q3, mediant, q1t, q3t;
         double trmean = 0.0, trmeant = 0.0;

         median = (double)(Nconv[r_in.N_data/2] + Nconv[(r_in.N_data-1)/2])/2.0;
        q1 = (Nconv[r_in.N_data/4] + Nconv[(r_in.N_data-1)/4])/2.0;
        q3 = (Nconv[3*r_in.N_data/4] + Nconv[(3*r_in.N_data-1)/4])/2.0;
               mediant = (double)(tconv[r_in.N_data/2] + tconv[(r_in.N_data-1)/2])/2.0;
        q1t = (tconv[r_in.N_data/4] + tconv[(r_in.N_data-1)/4])/2.0;
        q3t = (tconv[3*r_in.N_data/4] + tconv[(3*r_in.N_data-1)/4])/2.0;
        for(j=r_in.N_data/4; j<3*r_in.N_data/4; j++){
            trmean += Nconv[j];
            trmeant += tconv[j];
        }
        trmean /= (double)(3*r_in.N_data/4 - r_in.N_data/4);
        trmeant /= (double)(3*r_in.N_data/4 - r_in.N_data/4);
          
            printf("MC iterations to convergence: ");
            print_stats(stdout, n_conv_histogram->stats);
            printf("CPU time to convergence: ");
            print_stats(stdout, t_conv_histogram->stats);
          
            fprintf(stdout,
                    "\nN_data: %i, %i. \n <N_conv> = %f +- %f. \n <cpu_t_conv> = %f +- %f \n",
                    r_in.N_data, (int)n_conv_histogram->stats->N, n_conv_histogram->stats->mean, n_conv_histogram->stats->stddev_mean,
                    t_conv_histogram->stats->mean, t_conv_histogram->stats->stddev_mean);
            /* fprintf(fp_out3, */
            /*         "\nN_data: %i, %i. \n <N_conv> = %f +- %f. \n <cpu_t_conv> = %f +- %f .\n", */
            /*         r_in.N_data, (int)n_conv_histogram->stats->N, n_conv_histogram->stats->mean, n_conv_histogram->stats->stddev_mean, */
            /*         t_conv_histogram->stats->mean, t_conv_histogram->stats->stddev_mean); */
	 }
	 //}// end of loop over parameter to vary   
      
    for (j=0; j<r_in.N_data; j++){          
        permfree(&p1[j]); permfree(&p2[j]);
    }
    
    t_mainstop = GETCPUTIME;
    printf("number of path updates with log_p_a not(normal or zero): %i \n", n_path_p_a_nnoz);
    printf("number of lambda updates with log_p_a not(normal or zero): %i \n", n_lambdaIT_p_a_nnoz);
    printf("number of r/xi updates with log_p_a not(normal or zero): %i \n", n_rxi_p_a_nnoz);
    printf("Total CPU time (seconds): %f \n\n", t_mainstop - t_mainstart);
    printf("Clocks, calls, clocks/call: \n");
    printf("Run_around: %g %g %g \n", run_around_clocks, run_around_calls, run_around_clocks/run_around_calls);
    printf("Run_around_old: %g %g %g \n", run_around_old_clocks, run_around_old_calls, run_around_old_clocks/run_around_old_calls);
        //printf("get_cycle_start: %g %g %g \n", get_cycle_start_clocks, get_cycle_start_calls, get_cycle_start_clocks/get_cycle_start_calls); 
    printf("get_CD: %g %g %g \n", get_CD_clocks, get_CD_calls, get_CD_clocks/get_CD_calls);
    printf("get_CD_old: %g %g %g \n", get_CD_old_clocks, get_CD_old_calls, get_CD_old_clocks/get_CD_old_calls);
    printf("%g\n", (get_CD_clocks/get_CD_calls)/(get_CD_old_clocks/get_CD_old_calls));
    
    check_alloc_info(stdout);
   
    {
        int ijk; FILE* fpx;
        int cume=0;
        fpx = fopen("acc_etc.out", "w");

         fprintf(fpx, "# zeta_hist (dist of where trans are along path): \n");
        for(ijk=0; ijk<20; ijk++){
            fprintf(fpx, "%g %i %i \n", 0.05*(ijk+0.5), zeta_hist_a[ijk], zeta_hist_b[ijk]);
                //    fprintf(stdout, "%g %i %i \n", 0.05*(ijk+0.5), zeta_hist_a[ijk], zeta_hist_b[ijk]);
        }fprintf(fpx, "\n");
        
        fprintf(fpx, "# nstucks: \n");
        for(i=0; i<50; i++){
            fprintf(fpx, "%i  %i \n", i, nstucks[i]);
        }fprintf(fpx, "\n");
        
        fprintf(fpx, "# prop_lengths: \n");        
        for(ijk = 0; ijk<2*MAXPATHLENGTH; ijk++){
            cume += prop_lengths[ijk];
            fprintf(fpx, "%5i %5i %5i \n", ijk, prop_lengths[ijk], cume);
        }
        fprintf(fpx, "\n");

         fprintf(fpx, "# lpa_hist: \n");        
        for(ijk = 0; ijk<100; ijk++){
            fprintf(fpx, "%5i %g %i \n", ijk, 0.5*(ijk - 70 + 0.5), lpa_hist[ijk]);
        } fprintf(fpx, "\n");

        
   

        fprintf(fpx, "#track length i, ptl[i], atl[i]\n");
        {
            double sum_ptl = 0.0, sum_atl = 0.0;
            for(ijk=0; ijk<MAX_N_MARKERS_PER_CHROMOSOME; ijk++){
                sum_ptl += ptl[ijk]; sum_atl += atl[ijk];
            }
                for(ijk=0; ijk<40; ijk++){
                    fprintf(fpx, "%i %g %g \n", ijk, ptl[ijk]/sum_ptl, atl[ijk]/sum_atl);
                }
            }
     /*    printf("n_inv_prop: %i n_trans_prop: %i \n", n_inv_prop, n_trans_prop); */
/*           for(ijk=0; ijk<50; ijk++){ */
/*             fprintf(stdout, "%i %i \n", ijk, length_L_diff[ijk]); */
/*         } */
    }
    WAITFORKEY
        return 0;
}

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
    
