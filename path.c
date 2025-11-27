// path.c 

#include "grit.h"
#include "params.h"
#include "structs.h" 
#include "perm.h"
#include "perm_etc.h" 
#include "exts.h"
#include "path.h"
#include "path_etc.h"
#include "signs.h"
#include "rngs.h"
#include "gr_misc.h"


// #define NOTSTOPPROB (4.0*Epsilon*Epsilon)
#define NOTSTOPPROB (0.4*Epsilon) 
#define BETA 0.0 // Total prob of flat piece of q(l) for path update; BETA = 1 -> q(l) is flat, BETA = 0, no flat piece

#define PRODATON TRUE //Need to do prod_a_to_n for mcmcmc (doesn't cancel)
 
void append_step_to_path(Path* path, Step* step)  // append a step to a path 
{
    step->prev = path->end;
    path->end->next = step; 
    path->end = step; 
}  // end of function append_step_to_path
 

Reversal take_step_multi(Permutation* perm, Cycle_decomposition* the_cd, double Epsilon, double* step_prob, double lambda_ratio) 
{ 
        // takes a step & returns Reversal done, or {UNKNOWN, ...} if no step.
        //(i.e. if path has reached target perm and chose to stop.
        // value of *step_prob when this function called is irrelevant
    
    Reversal the_rev;

        //  printf("top of take_step_multi. \n");
  
    if((the_cd->n_mark - the_cd->n_chrom - the_cd->n_int_cycles + the_cd->c_g) == 0){ // reached permutation equivalent to target
        if(drand(&rng_long) >= NOTSTOPPROB){ // don't take another step
                //    printf("take_step_multi. stop branch\n");   
            the_rev = r_init;
            the_rev.ld = the_rev.rd = UNKNOWN;
            *step_prob = 1.0 - NOTSTOPPROB;
        }
        else{
                //     printf("take_step_multi. notstop branch. Epsilon, Notstopprob: %g   %g \n", Epsilon, NOTSTOPPROB);
                //  choose uar from all inversions and translocations
            the_rev = rand_reverse_r(perm, 1.0, step_prob, the_cd->elements);
            *step_prob *= NOTSTOPPROB; // step_prob was set to 1/(n_inv + n_trans) by rand_reverse_r on prev. line
            the_rev.delta_c_i = UNKNOWN; the_rev.delta_c_g = UNKNOWN;
        }
    }
    else{
        the_rev = delta_ci_cg_reverse_new(the_cd, Epsilon, perm->n_inv, perm->n_trans, step_prob, lambda_ratio);
        invert_translocate(perm, &the_rev, FALSE);
    }
    return the_rev; 
} // end of function take_step_multi

Path*
shorten_path(Path* path)
{
        // shortens the path by eliminating loops

    Step* s0, * s1, * s2, * s3;
    // Step* a1, * a2;
    Cycle_decomposition the_cd;
    Path* new_path;
    Step* new_step;
    double logprodAold, logprodAnew;
    
   
    #ifndef NDEBUG
     Step* a1 = path->start; Step* a2 = path->end;
    Permutation* p1 = copy_perm(a1->perm); Permutation* p2 = copy_perm(a2->perm);
    #endif
    
    
    
    printf("top of shorten path. path->length: %i \n", path->length); fflush(stdout);
    get_CD(path->start->perm, path->end->perm, &the_cd);
    printf("after initial get_CD in shorten..\n"); fflush(stdout);
    
    for(s0 = path->start, s1 = s0->next; s1 != path->end && s1 != NULL; s0 = s1, s1 = s1->next){
        for(s2 = s1->next; s2 != NULL; s2 = s2->next){
            if(are_perms_equivalent(s1->perm, s2->perm)){
                s0->next = s2; // go from s0 to s2, eliminating steps s1, s1->next, etc.
                while(s1 != s2){
                    s3 = s1;
                    s1 = s1->next;
                    path->length--;
                    path->length_a -= (1 + s3->n_dum); // 1 for the real inv/translocation, n_dum for the dummy events
                    if(s3->rev.is_inversion == TRUE) path->length_i--;
                    else path->length_t--;
                    path->length_f -= (s3->rev.is_fission == TRUE)? 1: 0;
                    permfree(&(s3->perm)); 
                    free(s3); n_step_alloc--;   
                }
            }
        }
    }
    assert(a1 == path->start); assert(a2 == path->end);
    assert(are_perms_equivalent(p1, path->start->perm)); assert(are_perms_equivalent(p2, path->end->perm)); 
    printf("in shorten... length after shortening: %i \n", path->length);
    printf("after main body of shorten_path.\n"); fflush(stdout);
    path->fission_factor = get_fission_factor(path);
    printf("after get_fission_factor in shorten_path.\n"); fflush(stdout);
    path->swap_end_factor = get_swap_end_factor(path);
    printf("after get_swap_end_factor in shorten_path.\n"); fflush(stdout);
    {
        int check_path_lens, check_path_consist;
        check_path_lens = (check_path_lengths(path) == TRUE);
        printf("in shorten... after check_path_lengths.\n"); fflush(stdout);
        check_path_consist = (check_path_consistency_multi(path) == 0); // checks that perms and reversals along path are consistent
        printf("in shorten...check path lengths, consist: %i %i \n", check_path_lens, check_path_consist); fflush(stdout);
    }

    new_step = step_alloc(copy_perm(path->start->perm)); // construct step whose perm is copy of path->start
    new_path = path_alloc(new_step);
    
    continue_path_with_new_distances(new_path, path->start, &logprodAold, &logprodAnew);
    check_path_next_prev(new_path);
    printf("in shorten_path, after continue... and check_path_next...\n");

    get_path_lengths(new_path);
    printf("in shorten. after continue... lengths: %i %i \n", path->length, new_path->length);
    freepath(path); // free the old path
    free(path); // 
    path = new_path;
    printf("path->length: %i \n", path->length);
    
   
          printf("checking shortened path: %i \n", check_path_selfconsistent(new_path));
              //  print_path(stdout, new_path);
      printf("bottom of shorten_path. \n"); fflush(stdout);
          //  getchar();
    return new_path;
    
}  // end shorten_path                
                
 
Path* gen_path(Permutation* p1, Permutation* p2, double Epsilon, double lambda_ratio, double limit_l, double* log_prob)
//               , double prob_dist) // generate a path from p1 to p2 
{
  Permutation* the_p;
 
    Step* the_step; 
    Cycle_decomposition the_cd; 
    Path* the_path;
    double Eps;
    int old_ci, old_cg;
    int dd; // delta_d, for checking
    int ninv = 0, ntrans = 0;
    double pathprob = 1.0;
    int dp1p2; // est dist from p1 to p2
    //double step_prob_check;
    
    *log_prob = 0.0; // initialize the path prob to 1
//    prob_dist = 1.0; // prob density for these distances = product(1/(l1*l2)), i.e. there is a factor for each step on path of
        // (l1*l2)^-1 where l1, l2 are the distances between the pairs of markers which are separated by the inversion/translocation
    
    the_p = copy_perm(p1);
        //   print_perm(stdout, the_p);  getchar();
    get_CD(the_p, p2, &the_cd);
    dp1p2 = the_cd.d;
        //     printf("d: %i \n", the_cd.d);
        // the 3rd argument in the following is the number of steps so far in the generated path 
    Eps = get_epsilon(&the_cd, Epsilon, 0, limit_l); // set epsilon using number of cycles to est. min length.

    the_step = step_alloc(the_p);
    the_path = path_alloc(the_step);
    the_path->n_stuck = 0;
    the_path->dd01_inv = FALSE; the_path->dd01_trans = FALSE;
    
    // #ifndef NDEBUG
    Permutation* old_p = the_p;
    // #endif
    the_p = copy_perm(the_p);
        //  printf("in gen_path, before loop\n");
    while ((the_step->rev = take_step_multi(the_p, &the_cd, Eps, &(the_step->prob), lambda_ratio)).is_inversion != UNKNOWN){
            //    print_rev(stdout, the_step->rev);
            //   print_perm(stdout, the_p);
            //     printf("in genpath. d: %i , the_step->rev.left %i \n", the_cd.d, the_step->rev.left);
        if(the_step->rev.A < 0.0){printf("in genpath. the_step->rev.A: %g \n", the_step->rev.A); getchar();}
          
        assert(check_rev_perm_consistent(the_step->perm, &(the_step->rev)) == TRUE);
        dd = the_step->rev.delta_c_g - the_step->rev.delta_c_i;
        if(the_step->rev.stuck == TRUE) the_path->n_stuck++;
        assert(real_close(step_prob_check = get_step_prob_multi(old_p, the_p, p2, &(the_step->rev), Eps, lambda_ratio), the_step->prob, DBL_EPSILON*RCFACTOR_1));
            //  printf("\n\n"); getchar();
            //  printf("step_probs: %g %g \n", the_step->prob, step_prob_check);
        pathprob *= the_step->prob;
        *log_prob += log(the_step->prob);
            //      printf("in genpath. step_prob, path_prob, log(path_prob), log_prob: %g %g %g %g \n", the_step->prob, pathprob, log(pathprob), *log_prob);
           
        the_path->length++;
        if(!the_step->rev.is_inversion){
            the_path->length_t++; // count steps which are translocations
            if(the_step->rev.is_fission) the_path->length_f++; // count steps which are fissions (translocations involving one initial empty chromosome)
        }
       
        if(rev_just_swaps_chromosome_ends(the_step->perm, the_step->rev) == FALSE){
            the_path->L++;
                //   printf("the_path->L: %i \n", the_path->L);
            if(!the_step->rev.is_inversion) { the_path->L_t++; } // count steps which are translocations
            else { the_path->L_i++; }
        }
        else{
            if(the_step->rev.is_fission) { printf("rev just swaps ends, and is fission\n");  }
        }
        the_path->prob *= the_step->prob;
        if(PRINT_DD_INV_TRANS && the_cd.d%5 == 0){
        /*     printf("d:%4i, steps:%4i is_inv: %i, dd: %i  n_ddm1: %4i %4i %4i,  n_dd0: %4i %4i %4i,  n_dd1: %4i %4i %4i\n", */
/*                    the_cd.d, the_path->length, the_step->rev.is_inversion, dd, */
/*                    the_cd.n_dd_inv[0], the_cd.n_dd_trans[0], the_cd.n_dd[0], */
/*                    the_cd.n_dd_inv[1], the_cd.n_dd_trans[1], the_cd.n_dd[1], */
/*                    the_cd.n_dd_inv[2], the_cd.n_dd_trans[2], the_cd.n_dd[2]); */
/*             printf("dci1 invs: %i  %i %i %i %i %i   dci1 trans:  %i  %i %i %i %i %i\n", */
/*                    the_cd.n_dci_inv[2], */
/*                    the_cd.n_dsw_inv_dci1[0], the_cd.n_dsw_inv_dci1[1], */
/*                    the_cd.n_dsw_inv_dci1[2], the_cd.n_dsw_inv_dci1[3], the_cd.n_dsw_inv_dci1[4], */
/*                    the_cd.n_dci_trans[2], */
/*                    the_cd.n_dsw_trans_dci1[0], the_cd.n_dsw_trans_dci1[1], */
/*                    the_cd.n_dsw_trans_dci1[2], the_cd.n_dsw_trans_dci1[3], the_cd.n_dsw_trans_dci1[4]); */
/*             printf("dcgm1 invs: %i  %i %i %i %i %i   dcgm1 trans:  %i  %i %i %i %i %i\n", */
/*                    the_cd.n_dcg_inv[0], */
/*                    the_cd.n_dsw_inv_dcgm1[0], the_cd.n_dsw_inv_dcgm1[1], */
/*                    the_cd.n_dsw_inv_dcgm1[2], the_cd.n_dsw_inv_dcgm1[3], the_cd.n_dsw_inv_dcgm1[4], */
/*                    the_cd.n_dcg_trans[0], */
/*                    the_cd.n_dsw_trans_dcgm1[0], the_cd.n_dsw_trans_dcgm1[1], */
/*                    the_cd.n_dsw_trans_dcgm1[2], the_cd.n_dsw_trans_dcgm1[3], the_cd.n_dsw_trans_dcgm1[4]); */
            printf("***d***: %i L, Li, Lt %i %i %i \n", the_cd.d, the_path->length, the_path->length - the_path->length_t, the_path->length_t);
             printf("ddm1 invs: %i  %i %i %i %i %i   ddm1 trans:   %i  %i %i %i %i %i\n",
                   the_cd.n_dcg_inv[0] + the_cd.n_dci_inv[2],
                    the_cd.n_dsw_inv_ddm1[0], the_cd.n_dsw_inv_ddm1[1],
                    the_cd.n_dsw_inv_ddm1[2], the_cd.n_dsw_inv_ddm1[3], the_cd.n_dsw_inv_ddm1[4],
                    the_cd.n_dcg_trans[0] + the_cd.n_dci_trans[2],
                    the_cd.n_dsw_trans_ddm1[0], the_cd.n_dsw_trans_ddm1[1],
                    the_cd.n_dsw_trans_ddm1[2], the_cd.n_dsw_trans_ddm1[3], the_cd.n_dsw_trans_ddm1[4]);

             printf("                    ddm1 nswbefore = 2 trans:   %i %i %i %i %i \n",
                    the_cd.n_ddm1_nsw2before[0], the_cd.n_ddm1_nsw2before[1], the_cd.n_ddm1_nsw2before[2],
                    the_cd.n_ddm1_nsw2before[3], the_cd.n_ddm1_nsw2before[4]);
               printf("             ddm1 nswbefore = 2 trans, ncomnochange:   %i %i %i %i %i \n",
                    the_cd.n_ddm1_ncomnochange[0], the_cd.n_ddm1_ncomnochange[1], the_cd.n_ddm1_ncomnochange[2],
                    the_cd.n_ddm1_ncomnochange[3], the_cd.n_ddm1_ncomnochange[4]); 
   
             printf("             ddm1 nswbefore = 2 trans, ncomdecrease:   %i %i %i %i %i \n",
                    the_cd.n_ddm1_ncomdecrease[0], the_cd.n_ddm1_ncomdecrease[1], the_cd.n_ddm1_ncomdecrease[2],
                    the_cd.n_ddm1_ncomdecrease[3], the_cd.n_ddm1_ncomdecrease[4]);
             printf("             ddm1 nswbefore = 2 trans, ncomdectozero:   %i %i %i %i %i \n",
                    the_cd.n_ddm1_ncomdectozero[0], the_cd.n_ddm1_ncomdectozero[1], the_cd.n_ddm1_ncomdectozero[2],
                    the_cd.n_ddm1_ncomdectozero[3], the_cd.n_ddm1_ncomdectozero[4]);
             
             printf("n_dsw_dcgm1: %i %i %i %i %i \n", the_cd.n_dsw_dcgm1[0],  the_cd.n_dsw_dcgm1[1],  the_cd.n_dsw_dcgm1[2],
                    the_cd.n_dsw_dcgm1[3],  the_cd.n_dsw_dcgm1[4]);
             printf("n_dsw_dci1: %i %i %i %i %i \n", the_cd.n_dsw_dci1[0],  the_cd.n_dsw_dci1[1],  the_cd.n_dsw_dci1[2],
                    the_cd.n_dsw_dci1[3],  the_cd.n_dsw_dci1[4]);
             printf("\n");
        }
        
            /*   printf("d, li, lt: %i %i %i \n", the_cd.d, the_path->length - the_path->length_t, the_path->length_t); */
/*         printf("n_dsw_inv_ddm1: %i %i %i %i %i \n", the_cd.n_dsw_inv_ddm1[0],  the_cd.n_dsw_inv_ddm1[1],  the_cd.n_dsw_inv_ddm1[2], */
/*                the_cd.n_dsw_inv_ddm1[3],  the_cd.n_dsw_inv_ddm1[4]); */
/*         printf("n_dsw_trans_ddm1: %i %i %i %i %i \n", the_cd.n_dsw_trans_ddm1[0],  the_cd.n_dsw_trans_ddm1[1],  the_cd.n_dsw_trans_ddm1[2], */
/*                the_cd.n_dsw_trans_ddm1[3],  the_cd.n_dsw_trans_ddm1[4]); */
/*         printf("n_dsw_dcgm1: %i %i %i %i %i \n", the_cd.n_dsw_dcgm1[0],  the_cd.n_dsw_dcgm1[1],  the_cd.n_dsw_dcgm1[2], */
/*                the_cd.n_dsw_dcgm1[3],  the_cd.n_dsw_dcgm1[4]); */
/*         printf("n_dsw_dci1: %i %i %i %i %i \n", the_cd.n_dsw_dci1[0],  the_cd.n_dsw_dci1[1],  the_cd.n_dsw_dci1[2], */
/*                the_cd.n_dsw_dci1[3],  the_cd.n_dsw_dci1[4]); */
/*         printf("\n");  */
        
        if(dd >=0){
            if(the_step->rev.is_inversion == TRUE) { ninv++; the_path->dd01_inv = TRUE; }
            else if(the_step->rev.is_inversion == FALSE) { ntrans++; the_path->dd01_trans = TRUE; }
        }
        old_ci = the_cd.n_int_cycles; old_cg = the_cd.c_g; // for checking
        get_CD(the_p, p2, &the_cd); // get post-step cd
            //     printf("d: %i \n", the_cd.d);
        assert(the_step->rev.delta_c_i == UNKNOWN || the_cd.n_int_cycles - old_ci == the_step->rev.delta_c_i); // check that ci, cg are what we specified
        assert(the_step->rev.delta_c_g == UNKNOWN || the_cd.c_g - old_cg == the_step->rev.delta_c_g);
     
        Eps = get_epsilon(&the_cd, Epsilon, the_path->length, limit_l); // set epsilon using number of cycles to est. min length.
       
        append_step_to_path(the_path, the_step = step_alloc(the_p)); //append new step
        old_p = the_p; // store pointer to this perm so still can access it after taking next step
        the_p = copy_perm(the_p); // copy it; copy will be transformed by next step
    } // end of main gen_path loop
        //   printf("in gen_path, after loop\n");
        //   printf("\n");
        //   print_perm(stdout,
        //   printf("in gen_path, after loop. \n");
        //  print_rev(stdout, the_step->rev);
   
    if(PRINT_DD_INV_TRANS) printf("d: %4i,  n_ddm1: %4i %4i %4i,  n_dd0: %4i %4i %4i,  n_dd1: %4i %4i %4i \n", the_cd.d, the_cd.n_dd_inv[0], the_cd.n_dd_trans[0], the_cd.n_dd[0],
               the_cd.n_dd_inv[1], the_cd.n_dd_trans[1], the_cd.n_dd[1],
               the_cd.n_dd_inv[2], the_cd.n_dd_trans[2], the_cd.n_dd[2]);

    the_path->prob *= the_step->prob; // last "step" has prob 1.0 - NOTSTOPPROB
    *log_prob += log(the_step->prob); // 
     
    assert(are_perms_equivalent_both(&the_cd, the_path->end->perm, p2)); // checks that have reached a perm (the_p) which is equivalent to the target (p2)
    permfree(&the_p); // last one, not used in path, so free it here.   
      
    the_path->length_a = the_path->length; // length_a will include dummies (see add_dummies), length only inversions & translocations
    the_path->length_i = the_path->length - the_path->length_t;

    nstucks[the_path->n_stuck]++;
        //   if(the_path->length != the_path->L)
        // printf("%i %i     %i %i    %i %i \n", the_path->length, the_path->L, the_path->length_i, the_path->L_i, the_path->length_t, the_path->L_t);
        //  printf("in genpath, check_LiLt: %i \n", check_LiLt(the_path));
    assert(check_LiLt(the_path) == TRUE);
    do_zeta_hist(the_path);

        //  printf("in gen_path. dp1p2, L, Li, Lt: %i  %i %i %i \n", dp1p2, the_path->L, the_path->L_i, the_path->L_t);
        //  getchar();
        //  printf("bottom of gen_path, \n");
        //   print_path(stdout, the_path);
    check_path_next_prev(the_path);
    return the_path;         
}  // end of function gen_path


void
add_dummies(Path* path, const Lambdas* L)
{
        // add dummy events to a path
        // each step has a lambdum = Lambda - (lambdaI*n_inv + lambdaT*n_trans);    
        // prob that step has n dummy events:  p(n) = a^n * (1 - a);
        // every step including path->end gets dummies
    
    Step* the_step = path->start, * the_next_step;
    double a; // dummy rate as fraction of total rate
    int done = FALSE;
 
        // loop over steps on path
    while(!done){
        done = (the_step == path->end);
        the_next_step = the_step->next;
	//	fprintf(stderr, "ZZZ: %ld %g   %ld %g \n", (long)the_step->perm->n_inv, L->lambdaI, (long)the_step->perm->n_trans, L->lambdaT) ;
            //   a = (L->Lambda - (the_step->perm->n_inv*L->lambdaI + the_step->perm->n_trans*L->lambdaT))/L->Lambda;
        a = lambdum_over_Lambda(the_step->perm, L); 
        the_step->n_dum = (int)(log(1.0 - drand(&rng_long))/log(a));
	//    printf("in add_dummies, a, n_dum: %g %i \n", a, the_step->n_dum);
         
        path->length_a += the_step->n_dum;  
        the_step = the_next_step; 
    }
} // end of function add_dummies


//inline
double
lambdum_over_Lambda(const Permutation* p, const Lambdas* L)
{
    return 1.0 - Lambda_real_over_Lambda(p, L);
}

double
log_prod_a_to_n(const Path* path, const Lambdas* L)
{
        // get log(prod(a^n)) = sum(n*log(a)), a = lambdum/Lambda; n = number of dummy events between real ones
        // sum is over steps on path including last step, i.e. path->end

    Step* step = path->start, * next_step;
    double result = 0.0;
    int done = FALSE;
    
        // loop over steps on path
    while(!done){
        done = (step == path->end);
        next_step = step->next;
        result += log(lambdum_over_Lambda(step->perm, L))*(double)step->n_dum;
        step = next_step; 
    }
    return result;
} // end of function log_prod_a_to_n

//inline
double
Lambda_real_over_Lambda(const Permutation* p, const Lambdas* L)
{
  //fprintf(stderr, "ZZZZZZZZZZ: %g %g %g %g %g\n",
  //	  p->A_i, L->lambdaI, p->A_t, L->lambdaT, L->Lambda);
    return (p->A_i*L->lambdaI + p->A_t*L->lambdaT)/L->Lambda;
}

double
log_prod_A(const Path* path)
{
        // go along a path;
        // at each step find A = l1*l2 (just looks at rev.A)
        // where l1, l2 are the lengths of the edges which get broken
        // return log(prod(A(i)));
    
    Step* step;
    double result = 1.0;
 
        // loop over steps on path
    for(step = path->start; step != path->end; step = step->next){
        result *= step->rev.A;
    }
    return log(result);
} // end of function log_prod_A

double
prod_1minusa(const Path* path, const Lambdas* L)
{
        // go along a path; get lambdums
        // get prod(Lambda_real/Lambda) over steps on path; Where Lambda_real is the total
        // rate of real events (i.e. inversions and translocations)
      
    Step* step = path->start, * next_step;
    double result = 1.0;
    int done = FALSE;
 
        // loop over steps on path
    while(!done){
        done = (step == path->end);
        next_step = step->next;
        result *= Lambda_real_over_Lambda(step->perm, L);
        step = next_step; 
    }
    return result;
} // end of function prod_1minusa


int
check_LiLt(const Path* path)
{
        // go along a path; find Li, Lt (the number of inversions and translocations
        // along path, not counting those that just swap chromosome ends
      
    Step* the_step = path->start, * the_next_step;
    int result;
    int done = FALSE;
    Cycle_decomposition the_cd;
    int Li=0, Lt=0;
    int just_swap_ends;

    if(the_step != path->end)
    {
            // loop over steps on path
        while(!done){
            the_next_step = the_step->next; 
            done = done || (the_next_step == path->end);

            just_swap_ends = rev_just_swaps_chromosome_ends(the_step->perm, the_step->rev);
            get_CD(the_step->perm, the_next_step->perm, &the_cd);
            if(the_cd.d > 0){
                if(the_step->rev.is_inversion) Li++;
                else Lt++;
                    //   assert(rev_just_swaps_chromosome_ends(the_step->rev) == FALSE);
                if(just_swap_ends != FALSE){
                    print_perm(stdout, the_step->perm);
                    print_rev(stdout, the_step->rev);
                    print_perm(stdout, the_next_step->perm);
                    printf("cd.d: %i  just_swap_ends: %i \n", the_cd.d, just_swap_ends); getchar();
                }
            }
            else{
                if(just_swap_ends != TRUE){
                    print_perm(stdout, the_step->perm);
                    print_rev(stdout, the_step->rev);
                    print_perm(stdout, the_next_step->perm);
                    printf("cd.d: %i  just_swap_ends: %i \n", the_cd.d, just_swap_ends);  getchar();
                } 
            }
            the_step = the_next_step; 
        }
    }
    result = (Li == path->L_i) && (Lt == path->L_t);
    if(result != TRUE){
            //  print_path(stdout, path);
        printf("%4i    %4i %4i %4i      %4i %4i %4i \n", result, path->length_i,  Li, path->L_i,   path->length_t, Lt, path->L_t);
        getchar();
    }
    return result;
} // end of function check_LiLt


void
get_LLiLt(Path* path)
{
        // go along a path; find Li, Lt (the number of inversions and translocations
        // along path, not counting those that just swap chromosome ends
        // just looks at revs
      
    Step* the_step = path->start, * the_next_step;
    int done = FALSE;

    path->L = 0;
    path->L_i = 0;
    path->L_t = 0;
    
    if(the_step != path->end)
    {
            // loop over steps on path
        while(!done){
            the_next_step = the_step->next; 
            done = done || (the_next_step == path->end);
            
            if(rev_just_swaps_chromosome_ends(the_step->perm, the_step->rev) == FALSE){
                path->L++;
                    //   printf("path->L: %i \n", path->L);
                if(!the_step->rev.is_inversion) { path->L_t++; } // count steps which are translocations
                else { path->L_i++; }
            }
          
            the_step = the_next_step; 
        }
    }
} // end of function get_LLiLt

void
do_zeta_hist(const Path* path)
{
        // quick and dirty
        // go along a path; for path with L steps, if step i (i=0,1,,L-1) is a translocation,
        // put zeta = (i+0.5)/L in histogram
        // to study whether translocations occur uniformly along generated paths
      
    Step* the_step = path->start, * the_next_step;
    int done = FALSE;
    int step_number = 0;

    if(the_step != path->end)
    {
            // loop over steps on path
        while(!done){
            the_next_step = the_step->next; 
            done = done || (the_next_step == path->end);

            if(the_step->rev.is_inversion == FALSE){
                    //  printf("L, step_number: %i %i \n", path->length, step_number);
                zeta_hist_a[(int)(20.0*((double)step_number + 0.5)/path->length)]++;
                trans_a++;
                if(!are_perms_equivalent(the_step->perm, the_next_step->perm)){
                    zeta_hist_b[(int)(20.0*((double)step_number + 0.5)/path->length)]++;
                    trans_b++;
                }
                n_trans_prop++;
            }
            else n_inv_prop++;
            step_number++;    
            the_step = the_next_step; 
        }
    }
        //  printf("zeta_hist_b[0] %i \n", zeta_hist_b[0]);
    if(step_number != path->length) {printf("step_number, path->length: %i %i \n", step_number, path->length); getchar();}
} // end of function zeta_hist
    
double
get_n_inv_avg(const Path* path)


{
     // go along a path, get avg(n_inv) over steps on path;
     // should last step be included? prob not, but this just affects efficiency of algorithm, not correctness
        // but if don't include last step then returns zero for zero length paths, which is bad
        // so include last step
  Step* the_step = path->start, * the_next_step;
  int n_steps = 0;
  double sum = 0.0;
  int done = FALSE;
 
      // loop over steps on path
  while (!done) {
      done =  (the_step == path->end);
      the_next_step = the_step->next;
      sum += (double)the_step->perm->n_inv;
      n_steps++;
      the_step = the_next_step; 
  }
   return sum/(double)n_steps;
} // end of get_n_inv_avg

double
get_A_i_avg(const Path* path)


{
     // go along a path, get avg(A_i) over steps on path;
     // should last step be included? prob not, but this just affects efficiency of algorithm, not correctness
        // but if don't include last step then returns zero for zero length paths, which is bad
        // so include last step
  Step* the_step = path->start, * the_next_step;
  int n_steps = 0;
  double sum = 0.0;
  int done = FALSE;
 
      // loop over steps on path
  while (!done) {
      done =  (the_step == path->end);
      the_next_step = the_step->next;
      sum += (double)the_step->perm->A_i;
      n_steps++;
      the_step = the_next_step; 
  }
   return sum/(double)n_steps;
} // end of get_A_i_avg

double
get_A_t_avg(const Path* path)


{
     // go along a path, get avg(A_i) over steps on path;
     // should last step be included? prob not, but this just affects efficiency of algorithm, not correctness
        // but if don't include last step then returns zero for zero length paths, which is bad
        // so include last step
  Step* the_step = path->start, * the_next_step;
  int n_steps = 0;
  double sum = 0.0;
  int done = FALSE;
 
      // loop over steps on path
  while (!done) {
      done =  (the_step == path->end);
      the_next_step = the_step->next;
      sum += (double)the_step->perm->A_t;
      n_steps++;
      the_step = the_next_step; 
  }
   return sum/(double)n_steps;
} // end of get_A_t_avg


double
get_n_trans_avg(const Path* path)
{
     // go along a path, get avg(n_trans) over steps on path;
        
  Step* the_step = path->start, * the_next_step;
  int n_steps = 0;
  double sum = 0.0;
  int done = FALSE;
 
      // loop over steps on path
  while (!done) {
       done =  (the_step == path->end);
       the_next_step = the_step->next;

       sum += (double)the_step->perm->n_trans;
      n_steps++;  
      the_step = the_next_step; 
  } 
  return sum/(double)n_steps;
} // end of get_n_trans_avg


double get_length_prob(const Run_info_in* r_in, int N, int L, int length) 
{ 
        // given a path of length L, returns prob. of choosing a 
        // sub_path of length "length" 
        // to propose a replacement for.
        // N in 1 chromosome case was number of markers - what should it be for M>1? perhaps N-M?
     
    static int first = TRUE; 
    static double pLl[MAXPATHLENGTH+1][MAXPATHLENGTH+1]; 
    static double Alpha = -100.0; 
    double probs[MAXPATHLENGTH+1], cume_probs[MAXPATHLENGTH+1];
    double the_prob, cume_prob, temp; 
    int j, l;

    if(L >= MAXPATHLENGTH){
        printf("In get_length_prob. Path length: %i \n", L);
    }
    /* if(L >= MAXPATHLENGTH - 10){ */
/*         printf("Path length: %i \n", L); */
/*         getchar(); */
/*     } */
    assert(L <= MAXPATHLENGTH);
    assert(length <= L); 
 
    if (UPDATE_WHOLE_PATH) { // just propose update to whole path each time! 
        return (L == length)? 1.0: 0.0;; // prob is 1.0 if L == length, zero otherwise 
    } 
     
    if (Alpha != r_in->Alpha){ // alpha value has been changed! (re)do first-time calculation 
        first = TRUE; 
        Alpha = r_in->Alpha;
    }
    
     
    if (first){
        if(r_in->L_prop_dist_type < 2 || r_in->L_prop_dist_type > 5){
            printf("In get_length_prob. Unrecognized L_prop_dist_type. Exitting \n");
            process_error();
        }
        cume_prob = 0.0; 
        for(j=0; j<=MAXPATHLENGTH; j++){
            if(r_in->L_prop_dist_type == 2){ // 1 - tanh shape
                double temp = r_in->Xi* ((double)j/((double)N*r_in->Alpha) - 1.0);
                the_prob = 0.5*(exp(-temp)/cosh(temp));
                    //  the_prob = 0.5*(1.0 - tanh(r_in->Xi* ((double)j/((double)N*r_in->Alpha) - 1.0) ) );
            }
            else if(r_in->L_prop_dist_type == 3){ // lorentzian shape
                temp = r_in->Xi* ((double)j/((double)N*r_in->Alpha) - 1.0);
                the_prob = 1.0/(1.0 + temp*temp);
            }
            else if(r_in->L_prop_dist_type == 4){ // lorentzian squared shape
                temp = r_in->Xi* ((double)j/((double)N*r_in->Alpha) - 1.0);
                the_prob = 1.0/((1.0 + 0.5*temp*temp)*(1.0 + 0.5*temp*temp));
            }
            else if(r_in->L_prop_dist_type == 5){ // gaussian shape
                temp = r_in->Xi* ((double)j/((double)N*r_in->Alpha) - 1.0);
                the_prob = exp(-(temp*temp));
            }
	    else{ // unknown L_prop_dist_type 
	      fprintf(stderr, "In get_length_prob; r_in->L_prop_dist_type has unrecognized value: %ld . Bye.\n",
		      (long)r_in->L_prop_dist_type);
	    }
            cume_prob += probs[j] = the_prob;
            cume_probs[j] = cume_prob;
        } 
 
        for(j=0; j<= MAXPATHLENGTH; j++){ 
            for (l=0; l<=j; l++){ 
                    // normalize, add flat piece
                pLl[j][l] = (1.0 - BETA)*probs[l]/cume_probs[j] + BETA/(double)(j+1);
            }
        } 
        first = FALSE;
    }
    return pLl[L][length]; 
}  // end of function get_length_prob
 
void get_rand_length(const Run_info_in* r_in, int N, int L, int* length, double* prob) 
{ 
        // get length of section to update, and prob of choosing that length 
        // L is length of main path
        // N is number of markers for now (should in depend on number of chromosomes  too?
    double rnumb, cume_prob, the_prob; 
    int l; 
 
    assert(L <= MAXPATHLENGTH); 
    if(L > MAXPATHLENGTH) { 
      fprintf(stdout, "In get_rand_length. path length: %i ; too long \n", L); 
      process_error();
    } 
    if(UPDATE_WHOLE_PATH){ 
        *length = L; 
        *prob = 1.0; 
    } 
     
    rnumb = drand(&rng_long); 
    cume_prob = 0.0; 
 
    for (l=0; l<=L; l++){ 
        the_prob = get_length_prob(r_in, N, L, l);
            //      printf("in get_rand_length: l, L, the_prob: %i %i %g \n", l, L, the_prob); 
        cume_prob += the_prob; 
        if (cume_prob > rnumb){ 
            *length = l; 
            *prob = the_prob; 
            return; 
        } 
    } 
    *length = L;
    *prob = get_length_prob(r_in, N, L, L); 
    return; 
} // end of function get_rand_length


void get_section(Path* the_path, int update_length, Path* subpath) 
{ 
// sets subpath->start, subpath->end to pointers to steps at the two ends of section to update.
// also sets subpath->start->prev to ptr to step before subpath->start 
   
    int i, start;
    Step* before_step = NULL;
 
        // start is the number of the first step in the section to be updated 
    start = (int)(drand(&rng_long)*(the_path->length - update_length + 1)); // int in range [0, the_path->length - update_length]
    subpath->start = the_path->start; 
       
    for (i=0; i<start; i++){
        before_step = subpath->start;
        subpath->start = subpath->start->next; 
    } 
    subpath->end = subpath->start; 
    for (i=0; i<update_length; i++){ 
        subpath->end = subpath->end->next; 
    }
    subpath->length = update_length;
#ifdef MYDEBUG
    if(before_step != subpath->start->prev)printf("before_step, subpath_start->prev: %i %i \n", before_step, subpath->start->prev);
#endif
    assert(before_step == subpath->start->prev);
} // end of function get_section


double get_path_prob_multi(const Path* the_path, const Permutation* targ_perm, const double Epsilon, double lambda_ratio, double* log_prob) 
{ 
        // steps along path, getting CD at each step 
        // and then getting prob. of taking each step 
        // returns path prob = product of all step probs, 
        // including a prob of (1.0 - NOTSTOPPROB) for 
        // deciding to stop

        // now (with end cycles and int cycles) the_path does not necessarily
        // end at targ_perm. CDs here should be rel to targ_perm
 
    Step* step, * next_step; 
    Cycle_decomposition cd1, cd2; 
    Cycle_decomposition* the_cd = &cd1, * next_cd = &cd2, * temp_cd; 
    double the_path_prob = 1.0, step_prob;
    int i, delta_c_int, delta_m, delta_c_g, delta_d;
    int is_same_eds = UNKNOWN;
    int waitatend = FALSE;
    int path_consistent;
    int nsw, nsw_next, dnsw;
    

        //   printf("top of get_path_prob_multi. path->length: %i \n", the_path->length);
        //   print_perm(stdout,the_path->start->perm); print_perm(stdout, the_path->end->perm); printf("***\n");
    *log_prob = 0.0; // initialze to path prob of 1.0
    step = the_path->start; 
 
        // is the following get_CD needed? - Yes
    get_CD(step->perm, targ_perm, the_cd);
    nsw = n_syn(step->perm, targ_perm);

    for (i=0; i < the_path->length; i++){ 
        next_step = step->next; 
        assert(next_step != NULL);
        
        get_CD(next_step->perm, targ_perm, next_cd);

        delta_c_g = next_cd->c_g - the_cd->c_g;
        delta_c_int = next_cd->n_int_cycles - the_cd->n_int_cycles;
        delta_m = next_step->perm->n_chrom_nonempty - step->perm->n_chrom_nonempty;
        delta_d = delta_c_g - delta_c_int;

         nsw_next = n_syn(next_step->perm, targ_perm);
         dnsw = nsw_next - nsw; 
        
        if((step->perm->n_chrom_nonempty == targ_perm->n_chrom_nonempty)
           && (the_cd->n_int_cycles == targ_perm->n_mark - targ_perm->n_chrom_nonempty)){ // step->perm equivalent to targ_perm
            step_prob = NOTSTOPPROB;
            step_prob /= (double)(step->perm->n_inv + step->perm->n_trans);
        }
        else{
            if(delta_d == 0){
                is_same_eds = find_whether_same_eds(step->perm, the_cd, &(step->rev));
                if(is_same_eds == -100){waitatend = TRUE;
                path_consistent = (check_path_consistency_multi(the_path) == 0); // checks that perms and reversals along path are consistent}
                }
            }
            assert(check_rev_perm_consistent(step->perm, &(step->rev)) == TRUE);
              
              /*   printf("in get_path_prob. perms, next perm, targ perm: \n"); */
/*                 print_perm(stdout, step->perm); */
/*                 print_perm(stdout, next_step->perm); */
/*                  print_perm(stdout, targ_perm); */
                     //   nsw1 = n_switches1(step->rev.marker_numbers, step->perm, targ_perm);
                     //   nsw2 = n_switches1(step->rev.marker_numbers, next_step->perm, targ_perm);
                    // nsw1 = n_syn(step->perm, targ_perm);
            
             /*    nsw_next = n_syn(next_step->perm, targ_perm); */
                
/*                 dnsw = nsw_next - nsw; */
                    //   printf("in get_path_prob. ndsw: %i \n", dnsw);
                    //      printf("in get_path_prob. dnsw: %i nsw1, nsw2: %i %i \n", dnsw, nsw1, nsw2);
                {
                  /*   int nsw1, nsw2; */
/*                     nsw1 =  n_syn(step->perm, targ_perm); */
/*                     nsw2 =  n_syn(next_step->perm, targ_perm); */
                        //  dnsw = nsw2 - nsw1;
                    step_prob = get_delta_d_prob_new(the_cd, delta_d, is_same_eds, step->rev.is_inversion,
                                                     step->perm->n_inv, step->perm->n_trans, Epsilon, lambda_ratio, dnsw);
                
                        //  printf("nsw, nsw_next, nsw1, nsw2: %i %i %i %i \n", nsw, nsw_next, nsw1, nsw2);
                }
            
        }
           
        assert(delta_m >= -1 && delta_m <= 1);
        assert(delta_c_int >= -1 && delta_c_int <= 1);
        assert(!(delta_m == -1 && delta_c_int == -1) && !(delta_m == 1 && delta_c_int == 1));
        
        the_path_prob *= step_prob;
        *log_prob += log(step_prob);
          /*      printf("in get_path_prob_multi. step_prob, the_path_prob, log(the_path_prob), *log_prob: %g %g %g %g \n", */
/*                       step_prob, the_path_prob, log(the_path_prob), *log_prob); */
       
            // swap which cd ptr points to which cd 
        temp_cd = the_cd; 
        the_cd = next_cd; 
        next_cd = temp_cd; 
       
        step = next_step;
        nsw = nsw_next;
            // printf("nsw: %i %i  %i %i \n", nsw, nsw_next, n_syn(step->perm, targ_perm), n_syn(next_step->perm, targ_perm));
    }
        //   printf("in get_path_prob... bottom1. path_prob, log(path_prob), log_prob: %g %g %g\n", the_path_prob, log(the_path_prob), *log_prob);
    the_path_prob *= 1.0 - NOTSTOPPROB; //probability of stopping upon reaching targ permutation - every path has this
        //  printf("in get_path_prob... bottom2. path_prob, log(path_prob), log_prob: %g %g %g\n", the_path_prob, log(the_path_prob), *log_prob);
    *log_prob += log(1.0 - NOTSTOPPROB); 
        //   printf("in get_path_prob... bottom3. path_prob, log(path_prob), log_prob: %g %g %g\n", the_path_prob, log(the_path_prob), *log_prob);

    if(waitatend == TRUE){
        printf("path_consistent: %i \n", path_consistent);
        print_path(stdout, the_path); getchar();
    }
    
    return the_path_prob; 
} // end of function get_path_prob_multi

double
rel_pi(State* state)
{
        // returns exp(-lambda)*(lambda/NpMchoose2)^L*(0.5^L_t)/L!;
        // i.e. the rel prob of path given lambda for case where lambdaI = 2*lambdaT and lambdas are fixed.
    Path* the_path = state->path;
    int i, L = the_path->length;
    int L_t = the_path->L_t;
    int n_edges = the_path->start->perm->n_mark + the_path->start->perm->n_chrom;
    double lambda = state->L->lambdaI*(double)(n_edges*(n_edges-1)/2); // 
    double result = exp(-lambda);
    double x = 2.0*lambda/(double)(n_edges*(n_edges-1));
       
    for(i=1; i<=L; i++){
        result *= x;
        result /= (double)i;
        if(i<= L_t) result *= 0.5;
    }
     if(!is_normal(result)){
        printf("in rel_pi result: %g, result is not normal.\n", result); 
    }
    return result;
} // end of rel_pi

double log_pi_ratio_indep(const Path* old_path, const Path* new_path, const int NpM, const Lambdas* L)
{
        // like log(pi_ratio_indep()), but uses logs thoughout to avoid overflow problems
    double result = 0.0;

    if(new_path->length > MAXPATHLENGTH){ printf("in log_pi_ratio_indep. new_path->length > MAXPATHLENGTH \n"); return 0.0; }
    
    result += (double)(new_path->length_i - old_path->length_i)*log(L->lambdaI/L->Lambda);
    result += log(L->lambdaT/L->Lambda)*(double)(new_path->length_t - old_path->length_t);
        //   if(PRODATON) result *= prod_a_to_n(new_path, L)/prod_a_to_n(old_path, L); // cancels with part of q_ratio - don't need to do
    if(!is_normal_or_zero(result)){
        printf("in log_pi_ratio_indep result: %g, result is not (normal or zero).\n", result);
    }
    return result;  
} // end of log_pi_ratio_indep


double log_pi_ratio_1(int l_old, int l_new, double lambda) 
{
        // gives log( (lambda^l_new/l_new!)/(lambda^l_old/l_old!) )
    int i; 
    double result = 0.0;
    double log_lambda = log(lambda);
    if (l_new > l_old){ 
        for(i=l_old+1; i<=l_new; i++){ result -= log((double)i); }
    } 
    else if(l_new < l_old) { // l_new < l_old
        for(i=l_new+1; i<=l_old; i++){ result += log((double)i); }
    }
    result += (double)(l_new - l_old)*log_lambda;
    if(!is_normal_or_zero(result)){ printf("in log_pi_ratio_1 result: %g, result is not (normal or zero).\n", result); }
    return result; 
}   // end of function log_pi_ratio_1


double log_pi_ratio(int l_old, int l_t_old, int l_new, int l_t_new, int NpM, double Lambda) 
{     
    int i; // NpM = perm->n_mark + perm->n_chrom;
    double result = 0.0;
    double log_lambda = log(Lambda/(double)(NpM*(NpM-1)/2));   
   
    if (l_new > l_old){ 
        for(i=l_old+1; i<=l_new; i++){ result += log_lambda - log((double)i); }
    } 
    else if(l_new < l_old) { // l_new < l_old
        for(i=l_new+1; i<=l_old; i++){ result -= log_lambda - log((double)i); }
    }
    
    result += (double)(l_t_old - l_t_new)*log(2.0); // this puts in a factor of 1/2 for each translocation
     
    assert(logs_real_close(result, log_pi_ratio_check(l_old, l_t_old, l_new, l_t_new, NpM, Lambda), DBL_EPSILON*RCFACTOR_1));
    if(!is_normal_or_zero(result)){ printf("in log_pi_ratio result: %g, result is not (normal or zero) .\n", result); }
    
    return result; 
}   // end of function log_pi_ratio

double
log_pi_ratio_check(int l_old, int l_t_old, int l_new, int l_t_new, int NpM, double Lambda)
{
    double NpM_choose_2 = (double)(NpM*(NpM-1)/2);
    return (log_factorial(l_old) - log_factorial(l_new)) + log(NpM_choose_2/Lambda)*(double)(l_old - l_new) + log(2.0)*(double)(l_t_old - l_t_new);
} // end of log_pi_ratio_check

// ******************

double
log_pi_ratio_gen(Path* p_old, Lambdas* L_old, Path* p_new, Lambdas* L_new)
{
        // like pi_ratio_gen, but works with logs to avoid overflows
        //   Path* p_old = old->path, * p_new = new->path;
        //   Lambdas* L_old = old->L, * L_new = new->L;
    double result = 0.0;
    double log_factorial_ratio = 0.0;
    int i;

    result += (L_old->Lambda - L_new->Lambda);
    result += log(L_new->Lambda)*p_new->length_a - log(L_old->Lambda)*p_old->length_a;
    result += log(L_new->lambdaI/L_new->Lambda)*p_new->length_i - log(L_old->lambdaI/L_old->Lambda)*p_old->length_i;
    result += log(L_new->lambdaT/L_new->Lambda)*p_new->length_t - log(L_old->lambdaT/L_old->Lambda)*p_old->length_t;
        //  result += log(prod_a_to_n(p_new, L_new)) - log(prod_a_to_n(p_old, L_old));
    result += log_prod_a_to_n(p_new, L_new) - log_prod_a_to_n(p_old, L_old);
    

    if(p_old->length_a < p_new->length_a){
        for(i = p_old->length_a+1; i <= p_new->length_a; i++){
            log_factorial_ratio += log((double)i);
        }
        result -= log_factorial_ratio;
    }
    else  if(p_old->length_a > p_new->length_a){
        for(i = p_new->length_a+1; i <= p_old->length_a; i++){
            log_factorial_ratio += log((double)i);
        }
        result += log_factorial_ratio;
    }
    else {}; // L_a's are equal; do nothing
    if(!is_normal_or_zero(result)){
        printf("in log_pi_ratio_gen result: %g, result is not (normal or zero).\n", result);
    }
    return result;
}// end log_pi_ratio_gen


double  
log_pi_ratio_gen1(int old_length, int old_length_a, int old_length_i, int old_length_t, double old_Lambda, double old_lambdaI, double old_lambdaT,  
                  int new_length, int new_length_a, int new_length_i, int new_length_t, double new_Lambda, double new_lambdaI, double new_lambdaT)
{
        // doesn't do the prod_a_to_n ratio
    double result = 0.0;  
    double log_factorial_ratio = 0.0;  
    int i;  

    result += (old_Lambda - new_Lambda);  
    result += log(new_Lambda)*new_length_a  - log(old_Lambda)*old_length_a;  
    result += log(new_lambdaI/new_Lambda)*new_length_i - log(old_lambdaI/old_Lambda)*old_length_i;  
    result += log(new_lambdaT/new_Lambda)*new_length_t - log(old_lambdaT/old_Lambda)*old_length_t;  
        //   result *= prod_a_to_n(p_new, L_new)/prod_a_to_n(p_old, L_old);  

    if(old_length_a < new_length_a){  
        for(i = old_length_a+1; i <= new_length_a; i++){  
            log_factorial_ratio += log((double)i);  
        }  
        result -= log_factorial_ratio;  
    }  
    else  if(old_length_a > new_length_a){  
        for(i = new_length_a+1; i <= old_length_a; i++){  
            log_factorial_ratio += log((double)i); 
        } 
        result += log_factorial_ratio; 
    } 
    else {}; // L_a's are equal; do nothing
    
        //    printf("in log_pi_ratio_gen1. result: %g   \n", result);
     if(!is_normal_or_zero(result)){
        printf("in log_pi_ratio_gen1 result: %g, result is not (normal or zero)\n", result);
    }
    return result; 
}// end log_pi_ratio_gen1 

// ***********************

 
int
get_path_lengths(Path* path) 
{ 
        // get path lengths (length, length_t, etc.) by counting steps until reach path->end
        // also set values of path->length, path->length_t, etc.
    int count = 0, t_count = 0, f_count = 0, d_count = 0; 
    Step* the_step;
   
    the_step = path->start; 
    while(the_step != path->end){
        count++;
        d_count += the_step->n_dum; // add # dummy events for this step
        if(!the_step->rev.is_inversion) t_count++;
        if(the_step->rev.is_fission) f_count++;
        the_step = the_step->next;
    }
    path->length = count;
    path->length_t = t_count;
    path->length_i = count - t_count;
    path->length_f = f_count;
    d_count += the_step->n_dum; // add in dummies events of last step (after last real inv/trans)
    path->length_a = path->length + d_count;
    return count; 
} // end of function get_path_lengths


// compare updating with and without updating to the end
int
update_path_2ways(State* state, const Run_info_in* r_in, double exponent)
{
    State copy_of_state = *state;
    Lambdas copy_of_L = *(state->L);
        unsigned long seed = rng_long;
    int x, pathssame;
        // printf("top of update_paths_2ways\n");
        //  printf("state->path, state->path->length: %i %i \n", state->path, state->path->length);
        //  print_path(stdout, state->path);
    copy_of_state.path = copy_path(state->path);
        //  printf("in update_path_2ways; after copy_path. compare_paths: %i \n", compare_paths(state->path, copy_of_state.path, &x));
    copy_of_state.L = &copy_of_L;
    
    seedrand(seed);
        //   printf("bef1 %g %g %g %g  %i\n", drand(&rng_long), drand(&rng_long), drand(&rng_long), drand(&rng_long), state->path->length);
    update_path(state, r_in, exponent, FALSE);
        //  printf("aft1 %g %g %g %g  %i\n", drand(&rng_long), drand(&rng_long), drand(&rng_long), drand(&rng_long), state->path->length);
    seedrand(seed);
        //  printf("bef2 %g %g %g %g  %i\n", drand(&rng_long), drand(&rng_long), drand(&rng_long), drand(&rng_long), copy_of_state.path->length);
    update_path(&copy_of_state, r_in, exponent, TRUE);
        //  printf("aft2 %g %g %g %g  %i\n", drand(&rng_long), drand(&rng_long), drand(&rng_long), drand(&rng_long), copy_of_state.path->length); 
    pathssame = compare_paths(state->path, copy_of_state.path, &x);
    if(pathssame != 0){
        printf("in update_path_2ways. compare_paths: %i \n", pathssame);
        getchar();
        print_perm(stdout, state->path->start->perm);
        getchar();
        print_perm(stdout, copy_of_state.path->start->perm);
        getchar();
        return 0;
    }
    return 1;
}

int update_path(State* state, const Run_info_in* r_in, double exponent, int update_to_end) //
{ 
        // chooses a section of the_path to update, proposes 
        // an alternative subpath, gets acceptance prob. 
        // if proposal accepted, replaces old subpath with new one 
        // returns 1 if proposal accepted, 0 if rejected. 

    Path* the_path = state->path;
    //Step* step_before_subpath_start; 
    int old_subpath_length; 
    int proposed_length, proposed_length_i, proposed_length_t, proposed_length_f, proposed_length_a;
    // int proposed_L, proposed_L_i, proposed_L_t;
    double old_subpath_length_prob;   
    double q_signs_ratio, l_prob_ratio;
    
    Path* new_subpath, old_subpath; // old_subpath is a Path rather than a Path* 
    Permutation* new_starting_perm; 
    int flips[NMARKERMAX+1], flips_by_number[NMARKERMAX+1];
    int NpM = the_path->start->perm->n_mark + the_path->start->perm->n_chrom;
    int N = the_path->start->perm->n_mark;
    double Epsilon;
    double log_pirat;
    double log_new_subpath_prob, log_old_subpath_prob;
    double log_q_ratio, log_p_accept;
    int accept;
    //double log_q_dist_old, log_q_dist_prop;
    int max_path_length = r_in->MC3_max_pathlength;
    //static int first = TRUE;
    //Path* new_down_path;
    //Path old_down_path;
    // Path* regen_path; 
    // int llll, cmprpaths;
    // unsigned long seed;
    int do_sign_flipping = ((r_in->Choose_signs > 0) && (drand(&rng_long) < FLIPSIGNSTEPPROB));

    update_to_end = update_to_end || r_in->Use_distances;

    if(drand(&rng_long) < PROB_DDM1_PREFINV_USE_DNSWITCHES){
        use_dnswitches = TRUE; ddm1_prefer_inversions = TRUE;
    }
    else{
        use_dnswitches = FALSE; ddm1_prefer_inversions = FALSE;
    }

    if(the_path->length > max_path_length) {
        printf("L, max_path_length, exponent: %i %i %g \n", the_path->length, max_path_length, exponent);
    }
    if(!USE2EPSILONS || drand(&rng_long) < BIGEPSILONFRACTION){
        Epsilon = r_in->Epsilon;
    }
    else{
        Epsilon = EPSILON_BIG_FACTOR*r_in->Epsilon;
    }
 
        // choose a section to propose an update for
        // old_subpath_length_prob is the prob that length of the section to be replaced is old_subpath_length
    
    get_rand_length(r_in, N, the_path->length, &old_subpath_length, &old_subpath_length_prob);
    get_section(the_path, old_subpath_length,/*  &step_before_subpath_start, */ &old_subpath);
    
     if(old_subpath_length != get_path_lengths(&old_subpath)){ process_error(); }

         // generate a new subpath to propose as a replacement 
    new_starting_perm = copy_perm(old_subpath.start->perm);
    
    if (do_sign_flipping){     //  in unsigned case, flip marker(s) of new_starting_perm
        q_signs_ratio = flip_signs(r_in, old_subpath.start->perm, new_starting_perm, old_subpath.end->perm, flips); 
            //  printf("q_signs_ratio: %g \n", q_signs_ratio);
    }
    else { q_signs_ratio = 1.0; }   // signed case, leave new_starting_perm as is 

     new_subpath = gen_path(new_starting_perm, old_subpath.end->perm, Epsilon, state->L->r, -1.0, &log_new_subpath_prob); // final arg < 0 -> get_epsilon doesn't adjust epsilon
     permfree(&new_starting_perm);
     if(state->Mode >= 2) add_dummies(new_subpath, state->L); 
     proposed_length = the_path->length + (new_subpath->length - old_subpath.length);
    prop_lengths[proposed_length]++; // 
    proposed_length_i = the_path->length_i + (new_subpath->length_i - old_subpath.length_i); 
    proposed_length_t = the_path->length_t + (new_subpath->length_t - old_subpath.length_t);
    proposed_length_f = the_path->length_f + (new_subpath->length_f - old_subpath.length_f);
    proposed_length_a = the_path->length_a + (new_subpath->length_a - old_subpath.length_a);
    
  /*   proposed_L = the_path->L + (new_subpath->L - old_subpath.L); */
/*     proposed_L_i = the_path->L_i + (new_subpath->L_i - old_subpath.L_i); */
/*     proposed_L_t = the_path->L_t + (new_subpath->L_t - old_subpath.L_t); */
    

        //     results are not identical if following line moved after checks of new_subpath. Why? Doesn't happen if turn compile optimization off
        //  old_subpath.prob = get_path_prob_multi(&old_subpath, new_subpath->end->perm, Epsilon); // prob of proposing old path given new

        // a few checks
    new_subpath->fission_factor = get_fission_factor(new_subpath);
    assert(check_path_lengths(new_subpath) == TRUE);
  
    assert((get_path_prob_multi(new_subpath, old_subpath.end->perm, Epsilon, state->L->r, &log_new_subpath_prob_check),
            logs_real_close(log_new_subpath_prob, log_new_subpath_prob_check, 1.0e-11)));
        // end of checks
    old_subpath.prob = get_path_prob_multi(&old_subpath, new_subpath->end->perm, Epsilon, state->L->r, &log_old_subpath_prob); // prob of proposing old path given new
  
    if(proposed_length >= max_path_length) {  // don't allow transitions to very long paths
        accept = FALSE;
    }   
    else{ // reasonable path length proposed, calculate q, pi ratios, etc.
        double new_subpath_length_prob;     
        
            // pi ratio
        if(state->Mode >= 2) {
            log_pirat = log_pi_ratio_1(the_path->length_a, proposed_length_a, state->L->Lambda);
            log_pirat += log_pi_ratio_indep(&old_subpath, new_subpath, NpM, state->L);;
            assert(logs_real_close(log_pirat, (log_pi_ratio_gen1(the_path->length, the_path->length_a, the_path->length_i, the_path->length_t,
                                                            state->L->Lambda, state->L->lambdaI, state->L->lambdaT,
                                                            proposed_length, proposed_length_a, proposed_length_i, proposed_length_t,
                                                            state->L->Lambda, state->L->lambdaI, state->L->lambdaT)), 8.0e-12) == TRUE);
            if(PRODATON){ log_pirat += log_prod_a_to_n(new_subpath, state->L) - log_prod_a_to_n(&old_subpath, state->L); }
        }
        else {
            log_pirat = log_pi_ratio(the_path->length, the_path->length_t, proposed_length, proposed_length_t, NpM, state->L->Lambda);
        }
            // end of pi ratio

         // proposal probability ratio:
        
            //    path_prob_ratio = old_subpath.prob/new_subpath->prob;  //    old/new
            //    q_ratio = l_prob_ratio*path_prob_ratio*q_signs_ratio;
             new_subpath_length_prob = get_length_prob(r_in, N, proposed_length, new_subpath->length);
            l_prob_ratio = new_subpath_length_prob/old_subpath_length_prob;
            log_q_ratio = log_old_subpath_prob - log_new_subpath_prob + log(l_prob_ratio) + log(q_signs_ratio);
           if(state->Mode >= 2) {
            if(PRODATON) log_q_ratio += log_prod_a_to_n(&old_subpath, state->L) - log_prod_a_to_n(new_subpath, state->L) ;
                 log_q_ratio += log(prod_1minusa(&old_subpath, state->L) / prod_1minusa(new_subpath, state->L));
		 // seed = rng_long;
                // *********** put in factor for the lengths **********
              if(/* r_in->Use_distances ||  */update_to_end){ 
                double log_q_dist_ratio, log_prod_A_old, log_prod_A_new;
                double logprodAolddown, logprodAnewdown;
                log_prod_A_old = log_prod_A(&old_subpath);
                log_prod_A_new = log_prod_A(new_subpath);
                log_q_dist_ratio = log_prod_A_old - log_prod_A_new;

                    // put in factors for the part of the path downstream of the updated section
                    // for it's new distances compared with old
                   
                continue_path_with_new_distances(new_subpath, old_subpath.end, &logprodAolddown, &logprodAnewdown);
                log_q_dist_ratio += logprodAolddown - logprodAnewdown;
                  
                    //   log_q_ratio += log_prod_A(&old_subpath) - log_prod_A(new_subpath);
                 if(r_in->Use_distances) log_q_ratio -= log_q_dist_ratio;
            }
              
        }
        if(!is_normal_or_zero(log_q_ratio)){
            printf("in update_path. log_q_ratio = %g . Not (normal or zero)\n", log_q_ratio);
            printf("%g %g %g %g \n", log_old_subpath_prob, log_new_subpath_prob, log(l_prob_ratio), log(q_signs_ratio));
            printf("%i %i %i %i  \n", the_path->length, old_subpath_length, proposed_length, new_subpath->length);
            printf("%g %g \n", old_subpath_length_prob, new_subpath_length_prob);
        }
            // end of proposal prob ratio

            // acceptance prob
	log_p_accept = log_pirat*exponent + log_q_ratio; // exponent is 1/T
          if(LTHEAT == TRUE){ //
                //    printf("in update_path, LTHEAT factor\n");
                //  log_p_accept -= the_path->length_t*LTHEATFACTOR*log(exponent);
         /*    double ltheatpiratio = f_Ltheat(1.0/exponent, proposed_length_t)/f_Ltheat(1.0/exponent, the_path->length_t); */
/*             printf("L_t, old, prop: %i %i  ltheatpiratio: %g , exponent: %g \n", the_path->length_t, proposed_length_t, ltheatpiratio, exponent); */
/*             printf("f_Ltheat, old, prop: %g %g \n",  f_Ltheat(1.0/exponent, the_path->length_t), f_Ltheat(1.0/exponent, proposed_length_t)); */
            if(OLDLTHEAT){
                log_p_accept += log(f_Ltheat(r_in->T_hot, 1.0/exponent, proposed_length_t)) - log(f_Ltheat(r_in->T_hot, 1.0/exponent, the_path->length_t));
            }
            else{
                log_p_accept += log(f_Ltheat1(1.0/exponent, LTHEATBETA, state->L->Lambda, state->L->lambdaT, proposed_length_t))
                    - log(f_Ltheat1(1.0/exponent, LTHEATBETA, state->L->Lambda, state->L->lambdaT, the_path->length_t));
            }
        }

            accept = ACCEPTALL || (log_p_accept > 0.0) || (log(drand(&rng_long)) < log_p_accept); 
            
        {
            int ilpa;
            ilpa = (int)(log_p_accept/0.5) + 70;
            ilpa = (ilpa < 0)? 0:ilpa;
            ilpa = (ilpa > 99)? 99:ilpa;
            if(ilpa >=0 && ilpa <100){
                lpa_hist[ilpa]++;
            }
        }
    }
    

    {
        int dL = proposed_length - the_path->length + 50;
        dL = (dL < 0)? 0: dL;
        dL = (dL > 99)? 99: dL;
        acc[dL][3]++; //
      
        if(new_subpath->dd01_inv == TRUE){dd01_inv_all++;}
        if(new_subpath->dd01_trans == TRUE){dd01_trans_all++;}
        if(new_subpath->dd01_inv == TRUE || new_subpath->dd01_trans == TRUE ){dd01_iort_all++;}
        if(new_subpath->dd01_inv == TRUE && new_subpath->dd01_trans == TRUE ){dd01_iandt_all++;}
        
        if(!is_normal_or_zero(log_p_accept)) n_path_p_a_nnoz++;

        if (accept){
                /* ************************************************************ */
                /* ******************  ACCEPT THE PROPOSAL!! ****************** */
                /* ************************************************************ */
            if(new_subpath->dd01_inv == TRUE){dd01_inv_acc++;}
            if(new_subpath->dd01_trans == TRUE){dd01_trans_acc++;}
            if((new_subpath->dd01_inv == TRUE) || (new_subpath->dd01_trans == TRUE)){dd01_iort_acc++;}
            if((new_subpath->dd01_inv == TRUE) && (new_subpath->dd01_trans == TRUE)){dd01_iandt_acc++;}
            
            if(log_p_accept >= 0.0) acc[dL][0]++;
            else acc[dL][1]++;

       //insert new_subpath in the_path in place of old_subpath
            the_path->accept = 1; // indicates proposal accepted
            {
                Step* step_before; //e = step_before_subpath_start;
                step_before = old_subpath.start->prev;
                if (step_before != NULL){ 
                    if(do_sign_flipping){ // flip signs as necessary upstream of updated subpath
                        get_flips_by_number(old_subpath.start->perm, flips, flips_by_number); 
                        fix_signs(the_path, old_subpath.start, flips_by_number); // flip signs on steps upstream of old_subpath.start
                    }           
                    assert(the_path->start != old_subpath.start);
                
                        // replace old_subpath with new_subpath ...
                    step_before->next = new_subpath->start;
                    new_subpath->start->prev = step_before;
                    
                } 
                else{ 
                    assert(the_path->start == old_subpath.start);
                        // replace old_subpath with new_subpath ...
                    the_path->start = new_subpath->start;
                }
            }

            if(update_to_end/*  || r_in->Use_distances */){ // replace everything from start of old_subpath all the way to end with new_subpath
                assert(the_path->end->next == NULL);
                new_subpath->end->next = the_path->end->next; // should be NULL
                new_subpath->end->rev = the_path->end->rev; // should be all UNKNOWN
                    // if(old_subpath.end != the_path->end) printf("old_subpath.end, the_path->end: %p %p \n", old_subpath.end, the_path->end);
                old_subpath.end = the_path->end; // now old_subpath goes all the way to endpoint genome
                the_path->end = new_subpath->end; // the_path now ends at version of endpoint genome at end of new_subpath
               
            }
            else{ // just replaces old_subpath with new_subpath
                new_subpath->end->next = old_subpath.end->next;
                new_subpath->end->rev = old_subpath.end->rev;
  
                if (new_subpath->end->next == NULL){       
                    assert(the_path->end == old_subpath.end); 
                    the_path->end = new_subpath->end;  
                } 
                else{
                    assert(the_path->end != old_subpath.end); 
                    old_subpath.end->next = NULL;
                    new_subpath->end->next->prev  = new_subpath->end; 
                }
            }

            the_path->length = proposed_length;
            the_path->length_i = proposed_length_i;
            the_path->length_t = proposed_length_t;
            the_path->length_f = proposed_length_f;
            the_path->length_a = proposed_length_a;
           
            get_LLiLt(the_path); // get correct values for L, L_i, L_t fields
            assert(check_LiLt(the_path));

                //   printf("in update_path. before get_fission_factor. \n");
            the_path->fission_factor = get_fission_factor(the_path);
            the_path->swap_end_factor = get_swap_end_factor(the_path);
              
              if(!check_path_lengths(the_path)) getchar();
            assert(check_path_lengths(the_path) == TRUE);
               
            freepath(&old_subpath); // free the memory of old path 
            free(new_subpath); n_path_alloc--; // don't free the steps on the path new_subpath, 
                //  but free the Path structure which has pointers to the start and end steps of new_subpath. 
                //  the steps of new_subpath have been incorporated into the_path
          
            assert(check_path_consistency_multi(the_path) == 0); // checks that perms and reversals along path are consistent
        } 
        else{  // reject the proposal - path is unchanged
            if(new_subpath->dd01_inv == TRUE){dd01_inv_rej++;}
            if(new_subpath->dd01_trans == TRUE){dd01_trans_rej++;}
            if((new_subpath->dd01_inv == TRUE) || (new_subpath->dd01_trans == TRUE)){dd01_iort_rej++;}
            if((new_subpath->dd01_inv == TRUE) && (new_subpath->dd01_trans == TRUE)){dd01_iandt_rej++;}
            
            acc[dL][2]++;
            the_path->accept = 0; // indicates rejection of proposal
            freepath(new_subpath); free(new_subpath); n_path_alloc--; // free mem allocated in gen_path
        }
    }
    if(the_path->length - the_path->L < 0) {printf("length < L: length: %i L: %i \n", the_path->length, the_path->L); getchar();}
    else length_L_diff[the_path->length - the_path->L]++;
      
    return the_path->accept; 
}      // end of function update_path
 
 
double 
get_epsilon(Cycle_decomposition* the_cd, double Epsilon, int path_length, double limit_l) 
{ 
        // if limit_l is positive, decreases epsilon when predicted path length gets too long 
        // the_cd is the cycle decomposition of, say p1 relative to p2. path_length is the path_length 
        // upstream of p1 plus downstream of p2. Uses cd to find the minimum total path length, and if 
        // this is too large (compared with limit_l), epsilon gets reduced. 
    int i; 
    double Lest; 
    double l_max = limit_l, sqrt_l_max = sqrt(l_max); 
    double result = Epsilon; 
     
    if(limit_l < 0.0){ 
        return result; 
    } 
    else{
            //   Lest = (double)(path_length + N_marker + 1 - the_cd->n_cycles);  // 1 chromosome version
        Lest = (double)(path_length + the_cd->n_mark - the_cd->n_chrom - the_cd->n_int_cycles + the_cd->c_g); // multi-chromosome version
        for(i=0; i<10; i++){ 
            if(Lest > l_max + i*sqrt_l_max) result *= 0.1; 
            else {return result;} 
        } 
        return result; 
    } 
}

int
compare_paths(Path* path1, Path* path2, int* count)
{
        // return -3 if paths are not between same (i.e. equivalent) endpoints
        // (this happens in unsigned case)
        // return 0 if paths are equivalent
        // 1 if path1 > path2
        // -1 if path1 < path2
        // -2 if can't determine which is greater
        // first look at lengths,
        // then (if lengths equal) look at #s of translocations,
        // then (if #s of translocation equal) look at permutations along path
        // ignores dummy events
    int temp;
    int starts_same, ends_same;
    Step* step1, * step2;

        //  printf("path1: \n"); print_path(stdout, path1); printf("path2: \n"); print_path(stdout, path2);
    *count = 0;
    starts_same = are_perms_equivalent(path1->start->perm, path2->start->perm);
    ends_same = are_perms_equivalent(path1->end->perm, path2->end->perm);
                                                                        
        // assert(starts_same);
    assert(ends_same);
    if(!(starts_same)){ printf("in compare paths, start perms are not same!\n");}
    if(!(ends_same)) { printf("in compare paths, end perms are not same!\n");} // can't compare paths
    if(!(starts_same && ends_same)) return -3;
    if(path1->length > path2->length) {
            //   printf("in compare_paths. length1,2: %i %i \n", path1->length, path2->length);
            //   print_path(stdout, path1); print_path(stdout, path2); getchar();
        return 1;}
    else if(path1->length < path2->length) {
            //   printf("in compare_paths. length1,2: %i %i \n", path1->length, path2->length);
            //   print_path(stdout, path1); print_path(stdout, path2); getchar();
        return -1;}
        // don't distinguish on basis of length_t 
  /*   else if(path1->length_t < path2->length_t) {//printf("length1,2: %i %i \n", path1->length, path2->length); */
/*         return 10;} */
/*     else if(path1->length_t > path2->length_t) {//printf("length1,2: %i %i \n", path1->length, path2->length); */
/*         return -10;} */
    else { // lengths equal, look at indiv perms
        assert(path1->length == path2->length);
            //   assert(path1->length_t == path2->length_t);
        step1 = path1->start; step2 = path2->start;
           
        while(step1 != path1->end){
            assert(are_perms_equivalent(step1->perm, step2->perm));
               
            if( rev_just_swaps_chromosome_ends(step1->perm, step1->rev) && rev_just_swaps_chromosome_ends(step2->perm, step2->rev)){
                temp = 0; // step1 and step2 are considered to be the same (both just relabel chromosome ends)
                    //   printf("xxx\n");
            }
            else{
                temp = compare_revs(step1->rev, step2->rev);
            }
            if(temp != 0) {
                    //   print_path(stdout, path1); print_path(stdout, path2); getchar();
                return temp;
            }
            step1 = step1->next;
            step2 = step2->next;
            (*count)++;
        }
        return 0; // paths are equal
    }
    printf("shouldn't get here, in compare_paths\n");
} // end of compare_paths


int compare_revs(Reversal rev1, Reversal rev2)
{
        // assumes rev1, rev2 are applied to equivalent permutations
        // 0 if "equal"
    int temp;

     if(!rev1.to_acbd) {
         assert(!rev1.is_inversion);
        temp = rev1.numbs[2]; rev1.numbs[2] = rev1.numbs[3]; rev1.numbs[3] = temp;
    }
    if(!rev2.to_acbd) {
        assert(!rev2.is_inversion);
        temp = rev2.numbs[2]; rev2.numbs[2] = rev2.numbs[3]; rev2.numbs[3] = temp;
    }

        // 12 34    21 43    34 12    43 21

    if( ( (rev2.numbs[0]==rev1.numbs[0]) && (rev2.numbs[1]==rev1.numbs[1]) && (rev2.numbs[2]==rev1.numbs[2]) && (rev2.numbs[3]==rev1.numbs[3]) )
       ||  ( (rev2.numbs[1]==rev1.numbs[0]) && (rev2.numbs[0]==rev1.numbs[1]) && (rev2.numbs[3]==rev1.numbs[2]) && (rev2.numbs[2]==rev1.numbs[3]) )
       ||  ( (rev2.numbs[2]==rev1.numbs[0]) && (rev2.numbs[3]==rev1.numbs[1]) && (rev2.numbs[0]==rev1.numbs[2]) && (rev2.numbs[1]==rev1.numbs[3]) )
       ||  ( (rev2.numbs[3]==rev1.numbs[0]) && (rev2.numbs[2]==rev1.numbs[1]) && (rev2.numbs[1]==rev1.numbs[2]) && (rev2.numbs[0]==rev1.numbs[3]) ) ){
        return 0; // revs are equivalent
    }
    else{
        return -2; // not equivalent
    }
} // end compare_revs


         
Path*
copy_path(Path* src_path)
{
  Path* targ_path = (Path*)chcalloc(1, sizeof(Path));
    int length, compare;
    Step* src_step = src_path->start;
    Step* targ_step, * targ_next_step;
    targ_step = (Step*)chcalloc(1, sizeof(Step)); n_step_alloc++;

    *targ_path = *src_path;
    
    targ_step->rev = src_step->rev;
    targ_step->prob = src_step->prob;
    targ_step->n_dum = src_step->n_dum; // copy the numbers of dummy events
    targ_step->perm = copy_perm(src_step->perm);
    targ_path->start = targ_step; 
    
    while(src_step->next != NULL){
       
      targ_next_step = (Step*)chcalloc(1, sizeof(Step)); n_step_alloc++;
        targ_next_step->rev = src_step->next->rev;
        targ_next_step->prob = src_step->next->prob;
        targ_next_step->n_dum = src_step->next->n_dum; // copy the numbers of dummy events
        targ_next_step->perm = copy_perm(src_step->next->perm);
        targ_step->next = targ_next_step;
       
        src_step = src_step->next;
        targ_step = targ_next_step;  
    }
    targ_step->next = src_step->next; // which will be NULL
    targ_path->end = targ_step;
    
    compare = compare_paths(src_path, targ_path, &length);
    if(compare != 0){
        printf("compare path & copy (0 <-> equivalency): %i \n", compare);
        printf("source path: \n");
        print_path(stdout, src_path);
        printf("copy of path: \n");
        print_path(stdout, targ_path);
    }

    { // run along path and set prev fields of step to point to right steps
        Step* step = targ_path->start, * next_step;

        step->prev = NULL;
        while((next_step = step->next) != NULL){
            next_step->prev = step;
            step = next_step;
        }
    }
        
    return targ_path;
} // end of copy_path

double
get_step_prob_multi(Permutation* p, Permutation* next_p, Permutation* targ_p, Reversal* rev, double Epsilon, double lambda_ratio)
{
    Cycle_decomposition the_cd, next_cd;
    int delta_c_int, delta_c_g; // delta_m,
    double step_prob;
    int delta_d, is_same_eds = UNKNOWN;

 /*    printf("* * * * * * * * * * * * * *\n"); */
/*     print_perm(stdout, p); */
/*     print_perm(stdout, next_p); */
/*     printf("* * * * * * * * * * * * * *\n"); */
        
    get_CD(p, targ_p, &the_cd);
    get_CD(next_p, targ_p, &next_cd);
             
    delta_c_int = next_cd.n_int_cycles - the_cd.n_int_cycles;
    // delta_m = next_p->n_chrom_nonempty - p->n_chrom_nonempty;
    delta_c_g = next_cd.c_g - the_cd.c_g;
    delta_d = delta_c_g - delta_c_int;
        //  printf("in get_step_prob_multi. delta c_g, c_i: %i %i \n", delta_c_g, delta_c_int);
        
    if((p->n_chrom_nonempty == targ_p->n_chrom_nonempty)
       && (the_cd.n_int_cycles == targ_p->n_mark - targ_p->n_chrom_nonempty)){ //p equiv to targ_p
        step_prob = NOTSTOPPROB;
        step_prob /= (double)(p->n_inv + p->n_trans);
            //     printf("  in get_step_prob_multi. notstopprob branch. NOTSTOPPROB, step_prob: %g %g \n", NOTSTOPPROB, step_prob);
    }
    else{
        if(delta_d == 0){
            is_same_eds = find_whether_same_eds(p, &the_cd, rev);
        }
            //     printf("in get_step_prob_multi. before get_delta_d_prob_new\n");
        assert(check_rev_perm_consistent(p, rev) == TRUE);

          {
                int nsw1, nsw2, dnsw;

/*                                 printf("in get_step_prob. perms, next perm, targ perm: \n"); */
/*                 print_perm(stdout, p); */
/*                 print_perm(stdout, next_p); */
/*                  print_perm(stdout, targ_p); */
/*                  print_rev(stdout, *rev); */
                     //   nsw1 = n_switches1(rev->marker_numbers, p, targ_p);
                     //   nsw2 = n_switches1(rev->marker_numbers, next_p, targ_p);

                   nsw1 = n_syn(p, targ_p);
                nsw2 = n_syn(next_p, targ_p);
                dnsw = nsw2 - nsw1;
                    //       printf("in get_step_prob. dnsw: %i nsw1, nsw2: %i %i \n", dnsw, nsw1, nsw2);
             
                step_prob = get_delta_d_prob_new(&the_cd, delta_d, is_same_eds, rev->is_inversion,
                                                 p->n_inv, p->n_trans, Epsilon, lambda_ratio, dnsw);
            }
              //  step_prob = get_delta_d_prob_new(&the_cd, delta_d, is_same_eds, rev->is_inversion, p->n_inv, p->n_trans,   Epsilon, lambda_ratio, rev);
    }
    return step_prob;
} // end get_step_prob_multi

/* double */
/* get_step_prob(Cycle_decomposition* the_cd, int delta_d, int is_same_eds, int is_inversion, int n_inv, int n_trans, double Epsilon, */
/*               double lambda_ratio, int dnsw) */
/* { */
/*     int nsw1, nsw2, dnsw; */
/*     nsw1 = n_switches(rev->marker_numbers, p, targ_p); */
/*     nsw2 = n_switches(rev->marker_numbers, next_p, targ_p); */
/*     dnsw = nsw2 = nsw1; */
            
/*     step_prob = get_delta_d_prob_new(&the_cd, delta_d, is_same_eds, rev->is_inversion, */
/*                                      p->n_inv, p->n_trans, Epsilon, lambda_ratio, dnsw); */
/*     return step_prob; */
/* }  */

int find_whether_same_eds(Permutation* perm, Cycle_decomposition* the_cd, Reversal* rev)
{
        // returns TRUE (FALSE) if rev is (is not) a same eds reversal with respect to the_cd
    int L1, L2, R1, R2; // L, R are positions of edges cut
    int L, R;
    Cycle_element* edge1, * edge2;

    L1 = (rev->numbs[0] != UNKNOWN)? perm->p[rev->numbs[0]].position/2: UNKNOWN;
    L2 = (rev->numbs[1] != UNKNOWN)? perm->p[rev->numbs[1]].position/2: UNKNOWN;
    R1 = (rev->numbs[2] != UNKNOWN)? perm->p[rev->numbs[2]].position/2: UNKNOWN;
    R2 = (rev->numbs[3] != UNKNOWN)? perm->p[rev->numbs[3]].position/2: UNKNOWN;
    if(L1 != UNKNOWN) L = L1;
    else if(L2 != UNKNOWN) L = L2;
    else  return FALSE; // edge is on empty chromosome; distinct eds's
    if(R1 != UNKNOWN) R = R1;
    else if(R2 != UNKNOWN) R = R2;
    else  return FALSE; // edge is on empty chromosome; distinct eds's
   
    if(((L1 != UNKNOWN && L2 != UNKNOWN) && (L1 != L2)) || ((R1 != UNKNOWN && R2 != UNKNOWN) && (R1 != R2))){
        printf("in find_whether_same_eds. \n");
        print_perm(stdout, perm);
        printf("numbs:: %i %i %i %i \n",rev->numbs[0],rev->numbs[1],
               rev->numbs[2],rev->numbs[3]);
    
        printf("L1, L2, R1, R2: %i %i %i %i \n", L1, L2, R1, R2);
        printf("about to ext from find_whether_same_eds. \n");
        getchar();
        return -100;
    }

        //   L = L1; R = R1;
    edge1 = the_cd->elements + L;
    edge2 = the_cd->elements + R;

    if(edge1->position != L) {printf("in find_whether_same_eds, inconsistency: L, edge1->position: %i %i \n",
                                     L, edge1->position); getchar();}
    if(edge2->position != R) {printf("in find_whether_same_eds, inconsistency: R, edge2->position: %i %i \n",
                                     R, edge2->position); getchar();}
  /*   printf("L,R, edge1,2,positions: %i %i %i %i . edge1,2 eds's: %i %i \n", */
/*                  L, R, edge1->position, edge2->position, edge1->eds, edge2->eds); */
    return (edge1->eds == edge2->eds);
} // end of find_whether_same_eds


void
track_lengths(const Path* path, double* the_ptl, double* the_atl)
{
        // go along path and at each step find the number of
        // inversions of each "track length" (i.e. number of markers inverted) which are possible
        // and accumulate in the_ptl (possible track lengths)
        // also accumulate in the_atl (actual track lengths) the track lengths of the inversions actually done
        // exclude translocations and inversion which invert all markers on a chromosome
    
    Step* the_step = path->start, * the_next_step;
    int done = FALSE;
    int n_chrom = path->start->perm->n_chrom;
    int j;
    int left, right;
    Marker_end* ends[4], * ordered_ends[4];
    int order[4];
    
    if(the_step != path->end)
    {
            // loop over steps on path
        while(!done){
            the_next_step = the_step->next; 
            done = done || (the_next_step == path->end);

            for(j=0; j<n_chrom; j++){
                possible_track_lengths_1chrom(the_step->perm->n_markers_on_chromosome[j], the_ptl);
            }
            if(the_step->rev.is_inversion && (rev_just_swaps_chromosome_ends(the_step->perm, the_step->rev) == FALSE)){
                (void)get_4ends(the_step->perm, the_step->rev.marker_numbers, ends, ordered_ends, order);
                left = ordered_ends[0]->position/2; right = ordered_ends[3]->position/2;
              /*   if(the_step->rev.right - the_step->rev.left + 1 == 0){ */
/*                     printf("left, right: %i %i \n", the_step->rev.left, the_step->rev.right); getchar(); */
/*                 } */
/*                 the_atl[the_step->rev.right - the_step->rev.left + 1] += 1.0; */
                the_atl[right - left] += 1.0;
                    //   printf("atl[0], atl[1]: %i %i %i %i  %i %i\n", atl[0], atl[1],the_atl[0], the_atl[1], the_step->rev.left, the_step->rev.right);
            }
            the_step = the_next_step; 
        }
    }
} // track_lengths


void
possible_track_lengths_1chrom(int chrom_size, double* the_ptl)
{
        // adds to the_ptl the numbers of inversions with each track length
        // possible for a chromosome with chrom_size markers
    int i;

    if(chrom_size > MAX_N_MARKERS_PER_CHROMOSOME){
        printf("in possible_track_lengths_1chrom. chrom_size > MAX_N_MARKERS_ON_CHROMOSOME. chrom_size: %i %i \n",
               chrom_size, MAX_N_MARKERS_PER_CHROMOSOME);
        getchar();
    }
    for(i=1; i< chrom_size; i++){
        the_ptl[i] += (double)(chrom_size + 1 - i);
    }
}

/* Path* */
/* gen_new_dist_path(const Step* olddownstart, const Step* newsubpathend, double* logprodAold, double* logprodAnew) */
/* { */
/*         // olddownstart is start of a path (the path downstream of subpath being updated), newstart has perm in same eq. class */
/*         // but with different distances. */
/*         // generate a path starting with newstart by applying the revs along the path starting with oldstart, */
/*         // but with new distances */
/*         // returns the new path */
    
/*     Permutation* next_p, * the_p = copy_perm(newsubpathend->perm); */
/*     const Step* olddownstep = olddownstart; Step * the_step; */
/*     Reversal rev ; */
/*     Path* the_path; // the new path being constructed */

/*     logprodAold = logprodAnew = 1.0;  */

/*     the_step = step_alloc(the_p);    */
/*     the_path = path_alloc(the_step); */
    
/*         //run along path */
/*     while(olddownstep->next != NULL){   */
/*         next_p = copy_perm(the_p); */
/*         rev = olddownstep->rev; */
/*         logprodAold += rev.A; */
/*         invert_translocate(next_p, &rev, TRUE); // does rev to next_p, with new distances */
/*         the_step->rev = rev; */
/*         logprodAnew += rev.A; */
/*         append_step_to_path(the_path, the_step = step_alloc(next_p)); */

/*         the_p = next_p; */
/*         olddownstep = olddownstep->next; // step along old downstream path */
/*     } */

/*     return the_path; */
/* } */


void
continue_path_with_new_distances(Path* newpath, Step* oldsubpathend, double* logprodAold, double* logprodAnew)
{
        // extend newpath downstream of updated section, where path should be same as old path in terms of marker order rearrangement,
        // differing only in distances
        // olddownstart is start of a path (the path downstream of subpath being updated), newstart has perm in same eq. class
        // but with different distances.
    
    Permutation* the_p = newpath->end->perm; // copy_perm(newpath->end->perm);
    Step* olddownstep = oldsubpathend;
    Step * the_step = newpath->end;
    Reversal rev;
    //int first = TRUE;
     
    *logprodAold = *logprodAnew = 0.0;

    if(the_step->prev != NULL){
        if(the_step != the_step->prev->next){
            printf("in continue... first step. step, step->prev->next: %i %i \n", the_step, the_step->prev->next);
            getchar();
        }
    } // getchar();

    assert(are_perms_equivalent(olddownstep->perm, the_p));
            
        //run along path
    while(olddownstep->next != NULL){  
        the_p = copy_perm(the_p);
        rev = olddownstep->rev;
        if(!check_rev(the_p->n_mark, &rev)) getchar();
        *logprodAold += log(rev.A);
    
        fix_rev1(the_p, &rev); // fix rev so it may be applied directly to the_p, i.e. find appropriate chrom ends,
            // and put the 4 marker ends in order
        invert_translocate(the_p, &rev, TRUE); // does rev to next_p, with new distances
        the_step->rev = rev;
            
        *logprodAnew += log(rev.A);
            //   printf("the_step: %p \n", the_step);
        append_step_to_path(newpath, the_step = step_alloc(the_p));
       
        if(the_step->prev->next != the_step) {
            printf("in continue...\n");
            printf("the_step->prev: %p  the_step->prev->next: %p  the_step: %p \n", the_step->prev, the_step->prev->next, the_step);
            getchar();
        }
        
        olddownstep = olddownstep->next; // step along old downstream path
        
            //    if(!are_perms_equivalent(olddownstep->perm, the_p)   )
     /*    { */
/*             print_rev(stdout, rev); */
/*             printf("old\n"); print_perm(stdout, olddownstep->perm); */
/*             printf("new\n"); print_perm(stdout, the_p); */
/*         }    */
        assert(are_perms_equivalent(olddownstep->perm, the_p));
        
        the_step->n_dum = olddownstep->n_dum;
     
    }
    assert(check_path_next_prev(newpath));
        //  printf("bottom of continue path ... \n");   
}

Path* regenerate_path(Path* old_path)
{
        // start with old_path->start, apply old_path->start->rev,
        // etc., regenerating the whole path, but with chromosome ends, position, orientation
        // consistent at each step in the regenerated in the regenerated path

    Step* old_step = old_path->start, * new_step;
    Permutation* the_perm = copy_perm(old_step->perm);
    Reversal rev = old_step->rev;
    Path* new_path;

    new_step = step_alloc(the_perm);
    new_path = path_alloc(new_step);
    printf("OLD perm \n"); print_perm(stdout, old_step->perm);
    printf("NEW perm \n"); print_perm(stdout, the_perm); 
    while(rev.is_inversion != UNKNOWN){ // i.e. this is not last step
        new_step->rev = rev;
        print_rev(stdout, rev);
        new_step->n_dum = old_step->n_dum;
        the_perm = copy_perm(the_perm);
        invert_translocate(the_perm, &rev, FALSE); // keep same distances!
        printf("check markend_dists2: %i \n", check_markend_dists2(the_perm));
        append_step_to_path(new_path, new_step = step_alloc(the_perm));
       
        old_step = old_step->next; // step along old path
         printf("OLD perm \n"); print_perm(stdout, old_step->perm);
         printf("NEW perm \n"); print_perm(stdout, the_perm);
        rev = old_step->rev;
    }
    printf("Bottom of regenerate path\n");
    return new_path; 
}

int
check_path_selfconsistent(Path* the_path)
{
        // checks that all perms on path are selfconsistent
        // presently doesn't check consistency of revs with perms
    Step* step = the_path->start;
    int result = TRUE, psc;
    int count = 0;

    while(step != NULL){
      psc = check_perm_selfconsistent(step->perm);
            printf("step: %i perm selfconsistent: %i \n", count, psc);
        step = step->next;
        if(!psc) result = FALSE;
        count++;
    }
    return result;
}

void
print_path_chroms_involved(const Path* path)
{
    Step* step = path->start;
    int i;
    while(step != NULL){
        printf("%i  ", step->rev.n_chroms_involved);
        for(i=0; i<step->rev.n_chroms_involved; i++){
            printf("%i ", step->rev.chroms_involved[i]);
        }printf("\n");
            
        step = step->next;
    }
}

void
update_path_order(Path* path)
{
    int count_trans = 0, trans_to_move = 1 + (int)(path->length_t*drand(&rng_long));
    Step* step, * start, * end;
    int back, forw;

    printf("in update_path_order. check_path_next_prev: %i \n", check_path_next_prev(path));
    step = path->start;
    if(!step->rev.is_inversion) count_trans++;
    while(count_trans < trans_to_move){
        step = step->next;
        if(!step->rev.is_inversion) count_trans++;
    }
        // now step points to a step with one of the translocations on path
    printf("length_t, trans_to_move: %i %i \n", path->length_t, trans_to_move);
    print_path_chroms_involved(path);
    commutation_path(path, step, &start, &end, &back, &forw);
    printf("back, forw: %i %i \n", back, forw);
}

void
commutation_path(Path* path, Step* step, Step** start, Step** end, int* backward, int* forward)
{
        // returns a subpath of path which contains step,
        // and which has property that step->rev commutes with
        // all the other revs on the subpath

    Step* other_step;

    *backward = 0; *forward = 0;

    printf("step chroms involved: %i %i \n", step->rev.chroms_involved[0], step->rev.chroms_involved[1]);
        // work backward
    other_step = step;
        //  printf("back step chroms involved: %i %i \n", other_step->rev.chroms_involved[0], other_step->rev.chroms_involved[1]);
    while(other_step->prev != NULL && revs_commute(&step->rev, &other_step->prev->rev)){
   
        other_step = other_step->prev;
           printf("back step chroms involved: %i %i \n", other_step->rev.chroms_involved[0], other_step->rev.chroms_involved[1]);  
           (*backward)++; printf("backward: %i \n", *backward);
           
           if(other_step->rev.chroms_involved[0] < 0) {
               print_path(stdout, path); print_perm(stdout, other_step->perm);
           }
           printf("after print_perm\n");
    }
    *start = other_step;

        // work forward
    other_step = step->next;
    printf("forw step chroms involved: %i %i \n", other_step->rev.chroms_involved[0], other_step->rev.chroms_involved[1]); 
    while(other_step->next != NULL && revs_commute(&step->rev, &other_step->rev)){
        other_step = other_step->next;
        printf("forw step chroms involved: %i %i \n", other_step->rev.chroms_involved[0], other_step->rev.chroms_involved[1]);
        (*forward)++; printf("forward: %i \n", *forward);
    }
    *end = other_step;
    printf("backward, forward: %i %i \n", *backward, *forward);
}

int revs_commute(Reversal* rev1, Reversal* rev2)
{
    int i,j;

    for(i=0; i<rev1->n_chroms_involved; i++){
        for(j=0; j<rev2->n_chroms_involved; j++){
            if(rev1->chroms_involved[i] == rev2->chroms_involved[j]) {
                printf("in revs_commute. chroms: %i %i \n", rev1->chroms_involved[i], rev2->chroms_involved[j]); return FALSE;
            }
        }
    }
    return TRUE;
}

int
check_path_next_prev(Path* path)
{
  //int s; //, spn, snp;
    Step* step = path->start;
    int result = TRUE;
    int count = 0;
    while(step->next != NULL){
      //s = (int)step;
        //snp = step->next == NULL? (int)NULL: (int)step->next->prev;
        //spn = step->prev == NULL? (int)NULL: (int)step->prev->next;
            //  printf("%i %i %i   %i\n", spn, step, snp,  (s = snp && s == spn)); 
        if(step->next->prev != step){
                   printf("%i  %i %i %i   %i\n", count, step->prev->next, step, step->next->prev,  (step == step->prev->next && step == step->next->prev)); 
            getchar();
            result = FALSE;
        }
        step = step->next; count++;
    }
    if(!result) { print_path_next_prev(path); getchar(); }
    return result;
}

void
print_path_next_prev(Path* path)
{
    Step* step = path->start;
   
    while(step->next != NULL){
        printf("step->prev, step, step->next:  %i %i %i    %i\n", step->prev, step, step->next, step->next->prev == step);
        step = step->next;
    }
}  

// end of path.c
