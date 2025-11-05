// run_chains.c

#include "grit.h"
#include "params.h"
#include "structs.h"
#include "perm.h"
#include "perm_etc.h" 
#include "path.h"
#include "path_etc.h"
#include "path_tree_etc.h"
#include "histogram.h"
#include "signs.h"
#include "chain.h"
#include "exts.h" 
#include "gamma_dev.h"
#include "rngs.h"
#include "gr_misc.h"
#include "run_chains.h"

int run_chains_new(Chain_set* the_chain_set, const Run_info_in* r_in, double* cpu_time)
{ 
        // updates the chains in the_chain_set many times
        // returns # of updates done to each chain

    int ii, iii, steps_so_far, runlength;
    double tstop, tstart, tnow; long tmin, tsec;
    Path_tree* path_tree;
    double p_accept_mc3;
    int n_mc3_swap[MAX_N_TEMPERATURES] = {0}, n_mc3_swap_try = 0;
    int mc3_max_pathlength;
    Chain* the_chain, * other_chain;
    State* the_state, * other_state; 
    Path* the_path;
    int last_print_n_steps = 20;
    double L_g = the_chain_set->chain[0][0]->state->path->start->perm->L_g;
    int converged = FALSE;
    double BoW_Burn_in = (r_in->Conv_srh - 1.0)*2.0;
    double BoW_L, BoW_Li, BoW_Lt, BoW_lambdaI, BoW_lambdaT;
    FILE* fp_paths = fopen("pathT0.out", "w");

    
        // ******************  initialization section ******************
    path_tree = alloc_path_tree();
    if(r_in->Burn_in > 0){ runlength = r_in->Burn_in*r_in->Stop_at; }
    else runlength = r_in->Max_runlength;
 
        //  ****************** end of initialization ******************   
        //  ****************** start iterating chains ******************
    
    tstart = GETCPUTIME;
    steps_so_far = 0;
      
    do{  // once through this loop does 1 update to each of N_temperatures*N_chain chains
        for (ii = 0; ii<r_in->N_chain; ii++){              
            for(iii=0; iii<r_in->N_temperatures; iii++){
                if(WRITERAW && (steps_so_far % WRITERAWEVERY == 0) && (ii == 0)) {
                    fprintf(fp_raw[iii], "\n %i ", steps_so_far);
                }
                the_chain = the_chain_set->chain[ii][iii];
                the_state = the_chain->state;
                the_path = the_state->path;
             
                update_chain(the_chain, r_in);
                if(iii==0 && (the_path->L <= 28 && the_path->L_t <= 11) && (steps_so_far % 100 == 0)) { print_path_brief(fp_paths, the_path); }
 
                    // write to Rawx.out files
                if(WRITERAW && (steps_so_far % WRITERAWEVERY == 0)){
                    if(iii == 0)track_lengths(the_path, ptl, atl); // add track lengths for this path to histograms ptl, atl
                    if(COUNT_WITHIN_EQUIVALENCE_CLASS_STEPS){
                        fprintf(fp_raw[iii], "%i %i %g %g ", (short)the_path->length_i, (short)the_path->length_t,
                                (float)the_state->L->lambdaI, (float)the_state->L->lambdaT);
                    }
                    else{
                        fprintf(fp_raw[iii], "%i %i %g %g ", (short)the_path->L_i, (short)the_path->L_t,
                                (float)the_state->L->lambdaI, (float)the_state->L->lambdaT);
                    }
                }

                if(iii == 0){ // keep track of number of hits for each short path
                    if(the_path->length <= MAXSHORTPATHLENGTH){   
                        path_tree_node_insert_multi(&(path_tree->start), the_state, 1);}
                }
            } // end of loop over temperatures (index iii)
        } // end of loop over N_chain chains (index ii)

             // do the metropolis-coupled mcmc swapping between chains at different temperatures
        if(r_in->Do_MC3_swap == TRUE) {
            for(ii=0; ii<r_in->N_chain; ii++){
                n_mc3_swap_try++;
                for(iii=0; iii<r_in->N_temperatures-1; iii++){
                    the_state = the_chain_set->chain[ii][iii]->state;
                          the_path = the_state->path; // needed? Yes.
                    other_state = the_chain_set->chain[ii][iii+1]->state;
                    the_chain = the_chain_set->chain[ii][iii];
                    other_chain = the_chain_set->chain[ii][iii+1];
                    p_accept_mc3 = exp(log_pi_ratio_gen(the_path, the_state->L, other_state->path, other_state->L)
                                       *(1.0/the_chain->temperature - 1.0/other_chain->temperature));
                    if(LTHEAT == TRUE){
                        if(OLDLTHEAT){
                            p_accept_mc3 /= f_Ltheat(r_in->T_hot, other_chain->temperature, other_state->path->length_t)
                                / f_Ltheat(r_in->T_hot, other_chain->temperature, the_path->length_t);
                            p_accept_mc3 /= f_Ltheat(r_in->T_hot, the_chain->temperature, the_path->length_t)
                                / f_Ltheat(r_in->T_hot, the_chain->temperature, other_state->path->length_t);
                        }
                        else{
                            p_accept_mc3 /= f_Ltheat1(other_chain->temperature, LTHEATBETA, other_state->L->Lambda, other_state->L->lambdaT, other_state->path->length_t)
                                / f_Ltheat1(other_chain->temperature, LTHEATBETA, the_state->L->Lambda, the_state->L->lambdaT, the_path->length_t);
                            p_accept_mc3 /= f_Ltheat1(the_chain->temperature,  LTHEATBETA, the_state->L->Lambda, the_state->L->lambdaT, the_path->length_t)
                                / f_Ltheat1(the_chain->temperature,  LTHEATBETA, other_state->L->Lambda, other_state->L->lambdaT, other_state->path->length_t);
                        }
                    }
                             
                    if(p_accept_mc3 > 1.0 || drand(&rng_long) < p_accept_mc3){ // swap between different temperature chains
                        the_chain_set->chain[ii][iii]->state = other_state;
                        the_chain_set->chain[ii][iii+1]->state = the_state;
                        n_mc3_swap[iii]++;
                    }
                }
            }
        } // end of if(r_in->Do_MC3_swap == TRUE) {
        
        chain_set_abw(the_chain_set);
            
        steps_so_far++;
        if(steps_so_far - last_print_n_steps >= 10 && steps_so_far > 1.02*(double)last_print_n_steps){ // print some stuff to stdout
            tnow = GETCPUTIME;
                //   printf("CPU time used: %g \n", tnow - tstart);
            tmin = (long)(tnow - tstart)/60;
            tsec =  (long)(tnow - tstart) - 60*tmin;
                //   printf("cpu time so far (min:sec): %i:%i \n", tmin, tsec);
            
            printf("\nUpdates: %i   CPU time so far (min:sec): %i:%i\n", steps_so_far, tmin, tsec);    
            if(r_in->N_temperatures > 1) {
                printf("    MCMCMC swaps: ");
                for(iii=0; iii<r_in->N_temperatures-1; iii++){
                    printf(" %i ", n_mc3_swap[iii]);
                } printf("Out of: %i", n_mc3_swap_try);
            }printf("\n");

            {
                double n_step = 0.0, n_acc_path = 0.0, n_acc_lambdaI = 0.0, n_acc_lambdaT = 0.0, n_acc_r = 0.0;
                printf("     T  P_a's:    path    lambdaI    lambdaT    r/xi \n"); 
                for(iii=0; iii<r_in->N_temperatures; iii++){
                    for(ii=0; ii<r_in->N_chain; ii++){
                        n_step = 0.0, n_acc_path = 0.0, n_acc_lambdaI = 0.0, n_acc_lambdaT = 0.0, n_acc_r = 0.0;
                        the_chain = the_chain_set->chain[ii][iii];
                        n_step += the_chain->N_steps_so_far;
                        n_acc_path += the_chain->n_acc_path;
                        n_acc_lambdaI += the_chain->n_acc_lambdaI;
                        n_acc_lambdaT += the_chain->n_acc_lambdaT;
                        n_acc_r += the_chain->n_acc_r;
                            //    printf("the_chain->n_acc_r, n_acc_r: %g  %g \n", the_chain->n_acc_r, n_acc_r);
                    }
                    the_chain_set->Pa_path[iii] = n_acc_path/n_step;
                    the_chain_set->Pa_lambdaI[iii] = n_acc_lambdaI/n_step;
                    the_chain_set->Pa_lambdaT[iii] = n_acc_lambdaT/n_step;
                    the_chain_set->Pa_r[iii] = n_acc_r/n_step;
                    printf("%8g     %9g  %9g %9g %9g \n", the_chain_set->chain[0][iii]->temperature, the_chain_set->Pa_path[iii],
                           the_chain_set->Pa_lambdaI[iii], the_chain_set->Pa_lambdaT[iii], the_chain_set->Pa_r[iii]);
                } // printf("\n");

            }
            
            printf("     T,   Avg:   L,        L_i,      L_t,      lambdaI,    lambdaT \n");
            for(iii=0; iii<r_in->N_temperatures; iii++){
                printf("%8g      %8g  %8g  %8g  %10g  %10g\n", the_chain_set->chain[0][iii]->temperature, 
                       the_chain_set->Avg_L[iii],
                       the_chain_set->Avg_Li[iii],
                       the_chain_set->Avg_Lt[iii],
                       the_chain_set->Avg_lambdaI[iii],
                       the_chain_set->Avg_lambdaT[iii]);
            }
               
            printf("     T,   B/W:   L,        L_i,      L_t,      lambdaI,    lambdaT \n"); 
            for(iii=0; iii<r_in->N_temperatures; iii++){
                printf("%8g      %8g  %8g  %8g  %10g  %10g\n", the_chain_set->chain[0][iii]->temperature, 
                       BoW_L = BoW(the_chain_set->B_L[iii], the_chain_set->W_L[iii]),
                       BoW_Li = BoW(the_chain_set->B_Li[iii], the_chain_set->W_Li[iii]),
                       BoW_Lt = BoW(the_chain_set->B_Lt[iii], the_chain_set->W_Lt[iii]),
                       BoW_lambdaI = BoW(the_chain_set->B_lambdaI[iii], the_chain_set->W_lambdaI[iii]),
                       BoW_lambdaT = BoW(the_chain_set->B_lambdaT[iii], the_chain_set->W_lambdaT[iii]));
                if(r_in->Burn_in < 0  &&  iii == 0 && !converged && steps_so_far >= r_in->Min_burn_in_length){
		  if( (converged = ((BoW_L < BoW_Burn_in) &&
                                    (BoW_Li < BoW_Burn_in) &&
                                    (BoW_Lt < BoW_Burn_in) &&
                                    (BoW_lambdaI < BoW_Burn_in) &&
                                    (BoW_lambdaT < BoW_Burn_in)) ) ){
                        runlength = r_in->Stop_at*steps_so_far;
                            //  runlength = (steps_so_far > r_in->Min_burn_in_length)? r_in->Stop_at*steps_so_far: r_in->Stop_at*r_in->Min_burn_in_length;
                        
                        
                        printf("*************  Burn-in done at %i updates.  ****************** \n", steps_so_far);
                            //  getchar();
                    }
                }
            }
                
          

            printf("    L_min: \n");
            for(iii=0; iii<r_in->N_temperatures; iii++){
                for(ii=0; ii<r_in->N_chain; ii++){ printf(" %3i ", the_chain_set->chain[ii][iii]->L_min);} printf("\n");
            } // printf("\n");
            printf("    L_t_min: \n");
            for(iii=0; iii<r_in->N_temperatures; iii++){
                for(ii=0; ii<r_in->N_chain; ii++){ printf(" %3i ", the_chain_set->chain[ii][iii]->L_t_min);} printf("\n");
            } // printf("\n");

          
            
            last_print_n_steps = steps_so_far;
        }
    }
    while (steps_so_far < runlength && !sigint_raised); // if sigint_raised (i.e. if control-c has been pressed), bail out here
    tstop = GETCPUTIME; 
    *cpu_time = tstop - tstart;
        // ****************** end of chain iteration section ******************
    
/*   fprintf(fp_out4, "Paths with length <= %i, and number of hits for each: \n", MAXSHORTPATHLENGTH); */
/*         print_path_tree(fp_out4, path_tree, TRUE); */
 
    return steps_so_far;
} // end of run_chains_new

double BoW(double B, double W)
{
    if(B == 0.0){
        return 0.0;
    }
    else return B/W;
}


int
fixed_lambda_mode(int mode)
{
    if(mode == 0 || mode == 2 || mode == 4) return TRUE;
    else return FALSE;
}

Histogram*
make_scaled_histogram(double xmin, double xrange, double minbw, double* binwidth)
{
#define MAXBINS 230
    double bw = 1.0;
    double nbins = xrange;

    while(bw > 10.0*minbw && nbins < MAXBINS){ // make nbins bigger than MAXBINS
        nbins *= 10.0;
        bw /= 10.0;
    }   

    if(xmin <0.0) xmin = 0.0;
    while(nbins > MAXBINS){
        nbins /= 2.0; bw *= 2.0;
        if (nbins > MAXBINS){
            nbins /= 1.25; bw *= 1.25;
            if (nbins > MAXBINS){
                nbins /= 2.0; bw *= 2.0;
            }
        }
    }
    printf("nbins, xmin, bw: %g %g %g \n", nbins, xmin, bw);
    *binwidth = bw;
    return create_histogram((int)nbins+4, bw*((int)(xmin/bw) - 0.5), bw);

}

int
get_histogram_parameters(double xmin, double nbins, int Integers, double* xmin_out, double* bw)
{
#define MAXBINS 230
    double minbw = 0.0;
        // minbw = (Integers == TRUE)? 1.0: 0.0;
    
    *bw = 1.0;
    if(Integers == TRUE){
         *xmin_out = (*bw)*((int)(xmin/(*bw)) - 0.5); // bins are centered on integers
    }
    else{
        while(*bw > 10.0*minbw && nbins < MAXBINS){ // make nbins bigger than MAXBINS to start
            nbins *= 10.0;
            *bw /= 10.0;
        }   

        if(xmin <0.0) xmin = 0.0;
        while(nbins > MAXBINS){
            nbins /= 2.0; *bw *= 2.0;
            if (nbins > MAXBINS){
                nbins /= 1.25; *bw *= 1.25;
                if (nbins > MAXBINS){
                    nbins /= 2.0; *bw *= 2.0;
                }
            }
        }
        *xmin_out = (*bw)*((int)(xmin/(*bw)));  // bin edges are multiples of bw
    }
    printf("nbins, xmin, bw: %g %g %g \n", nbins, xmin, *bw);
    return (int)nbins + 4;

}
        
