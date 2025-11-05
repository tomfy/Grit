// chain.c

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

//#include "math.h" <math.h> is included in grit.h
//#include <math.h>

#define NODISTLL (n_edges*(n_edges-1))

double srh(double a)
{
    return sqrt(1.0 + a);
}

void initialize_lambdas(State* state, const Run_info_in* r_in)
{
        // get values for lambdas for first element of chain
        //
    Path* path = state->path;
    Permutation* p1 = path->start->perm;
    int NpM = path->start->perm->n_edge;
    double L_g = state->path->start->perm->L_g; // genome length (== N+M for no distance case)
    int NpMchoose2 = (NpM-1)*NpM/2;
    double Iavg = get_n_inv_avg(path), Tavg = get_n_trans_avg(path), Ai_avg = get_A_i_avg(path), At_avg = get_A_t_avg(path);
        //   printf("Iavg, Tavg: %g %g \n", Iavg, Tavg);

    if(r_in->Use_distances == FALSE){
        L_g = (double)p1->n_edge;
        printf("L_g %g \n", L_g);
    }
     printf("L_g %g \n", L_g);

     // printf("XXXXX: in init lambdas. Update_lambdas_mode: %ld \n", (long)r_in->Update_lambdas_mode); getchar();
        //  printf("r_in->lambda_max, Lambda_max: %g %g \n", r_in->lambda_max, r_in->Lambda_max);
        
        // Modes 0, 1, no uniformization, no dummies added. (this calculation assumes r=2)
    if(r_in->Update_lambdas_mode == 0){ // fixed Lambda
        state->L->Lambda = r_in->Lambda;
        state->L->lambdaI = state->L->Lambda/(double)NpMchoose2;
        state->L->lambdaT = state->L->lambdaI/2.0;
    }
    else if(r_in->Update_lambdas_mode == 1){ // simulate Lambda
        state->L->Lambda = gamma_dev(path->length + 1, r_in->lambda_max*(double)NpMchoose2); //
        state->L->lambdaI = state->L->Lambda/(double)NpMchoose2;
        state->L->lambdaT = state->L->lambdaI/2.0;
    }
        // Modes 2, 3, use uniformization, dummies added. lambdaI, lambdaT parametrization
    else if(r_in->Update_lambdas_mode == 2){ // lambdaI, lambdaT fixed, with lambdaI/lambdaT = 2
        state->L->lambdaI = r_in->Lambda/(double)NpMchoose2;
        state->L->lambdaT = state->L->lambdaI/2.0;
        state->L->Lambda = get_Lambda(L_g, state->L->lambdaI, state->L->lambdaT);
        state->L->xi = get_xi(state->L->lambdaI, state->L->lambdaT); 
    }
    else if(r_in->Update_lambdas_mode == 3/*  || r_in->Update_lambdas_mode == 7 */){ // simulate from path to get init values of lambdaI
        state->L->lambdaI = gamma_dev(path->length_i + 1, r_in->lambda_max*Ai_avg)/Ai_avg;
        state->L->lambdaT = gamma_dev(path->length_t + 1, r_in->lambda_max*At_avg)/At_avg;
        state->L->Lambda = get_Lambda(L_g, state->L->lambdaI, state->L->lambdaT);
        state->L->xi = get_xi(state->L->lambdaI, state->L->lambdaT);
    }
     // Modes 4,5,6, use uniformization, dummies added. Lambda, r parametrization
    else if(r_in->Update_lambdas_mode == 4){ // Lambda fixed, r=2 fixed
        state->L->lambdaI = r_in->Lambda/(double)NpMchoose2;
        state->L->lambdaT = state->L->lambdaI/2.0;
        state->L->Lambda = get_Lambda(L_g, state->L->lambdaI, state->L->lambdaT);
        state->L->xi = get_xi(state->L->lambdaI, state->L->lambdaT); 
    }
     else if(r_in->Update_lambdas_mode == 5){ // r = 2 fixed, Lambda simulated
        state->L->lambdaI = gamma_dev(path->length_i + 1, r_in->lambda_max*Ai_avg)/Ai_avg;
        state->L->lambdaT = state->L->lambdaI/2.0;
        state->L->Lambda = get_Lambda(L_g, state->L->lambdaI, state->L->lambdaT);
       state->L->xi = get_xi(state->L->lambdaI, state->L->lambdaT);  
    }
    else if(r_in->Update_lambdas_mode >= 6){ // Lambda, r simulated
             //    printf("in init lambdas. length_i,t: %i %i \n", path->length_i, path->length_t);
             //    printf("lambda_max, Iavg: %g %g \n", r_in->lambda_max, Iavg);
             //   state->L->lambdaI = gamma_dev(path->length_i + 1, r_in->lambda_max*Iavg)/Iavg;
        
              state->L->lambdaI = (double)(path->length_i + 1)/Iavg;
              printf("lambdaI init no dist: %g \n", state->L->lambdaI);
               state->L->lambdaI = (double)(path->length_i + 1)/path->start->perm->A_i;
               printf("lambdaI init dist: %g \n", state->L->lambdaI);
               
               
               if(state->L->lambdaI > r_in->lambda_max) state->L->lambdaI = 0.9*r_in->lambda_max;
               if(Tavg > 0){
                       //   state->L->lambdaT = gamma_dev(path->length_t + 1, r_in->lambda_max*Tavg)/Tavg;

                   state->L->lambdaT = (double)(path->length_t + 1)/Tavg;
                   printf("lambdat init no dist: %g \n", state->L->lambdaT);
                   state->L->lambdaT = (double)(path->length_t + 1)/path->start->perm->A_t;
                   printf("lambdat init dist: %g \n", state->L->lambdaT);
                   printf("ITavg, Ai, At: %g %g %g %g \n", Iavg, Tavg, path->start->perm->A_i, path->start->perm->A_t);
                       //  getchar();
                   
                   if(state->L->lambdaT > r_in->lambda_max) state->L->lambdaT = 0.9*r_in->lambda_max;
               }
               else{ // 1 chromosome, no translocations
                   state->L->lambdaT = 0.0;
               }
              
                   //  printf("lambdaIT: %g %g \n", state->L->lambdaI, state->L->lambdaT);
             //  state->L->lambdaI = (double)path->length_i;
             //  state->L->lambdaT = (double)path->length_t;
               printf("L_g, lambdaI, lambdaT: %g %g %g    ", L_g, state->L->lambdaI, state->L->lambdaT);
        state->L->Lambda = get_Lambda(L_g, state->L->lambdaI, state->L->lambdaT);
        state->L->xi = get_xi(state->L->lambdaI, state->L->lambdaT);
        printf("Li, Lt: %i %i  Iavg, Tavg: %g %g \n", state->path->length_i, state->path->length_t, Iavg, Tavg);
             printf("in initialize_lambdas. Lambda, lambdaI, lambdaT, r, xi: %g %g %g %g %g\n",
                     state->L->Lambda ,state->L->lambdaI,state->L->lambdaT, state->L->lambdaI/state->L->lambdaT, state->L->xi);
    }else{
      fprintf(stderr, "Update_lambdas_mode has unknown value: %ld. Bye.\n",
	      (long)r_in->Update_lambdas_mode); exit(1);
    }
    state->L->r = state->L->lambdaI/state->L->lambdaT;
        //  printf("bottom of initialize lambdas\n");
} // end of initialize_lambdas



double
get_Lambda_old(int n_edges, double lambdaI, double lambdaT)
{
        // return a rate Lambda which is upper bound on rate of inversions plus translocations given LambdaI, LambdaT
        // so that dummy rate lambdum = Lambda - (lambdaI*n_inv + lambdaT*n_trans) is guaranteed to be positive.
    double alpha = LAMBDA_DEF_FACTOR*(double)NODISTLL;
    if(lambdaT > lambdaI){ 
        alpha *= lambdaT;
    }
    else { 
        alpha *= lambdaI;
    }
    return alpha;
} // end of get_Lambda

double
get_Lambda(double L_g, double lambdaI, double lambdaT)
{
        // return a rate Lambda which is upper bound on rate of inversions plus translocations given LambdaI, LambdaT
        // so that dummy rate lambdum = Lambda - (lambdaI*n_inv + lambdaT*n_trans) is guaranteed to be positive.
    double alpha = LAMBDA_DEF_FACTOR*L_g*L_g;
    if(lambdaT > lambdaI){ 
        alpha *= lambdaT;
    }
    else { 
        alpha *= lambdaI;
    }
    return alpha;
} // end of get_Lambda

double
get_xi(double lambdaI, double lambdaT)
{
        // return a value of xi corresponding to (lambdaI, lambdaT)
   
    if(lambdaT >= lambdaI){ 
      return lambdaI/lambdaT; // <= 1
    }
    else { 
      return 2.0 - lambdaT/lambdaI; // between 1 and 2
    }
} // end of get_xi

void
get_lambdaIT_from_Lambda_r_old(int n_edges, double Lambda, double r, double* lambdaI, double* lambdaT)
{
        // inverse of get_Lambda_alt
    double alpha = LAMBDA_DEF_FACTOR*(double)NODISTLL;
    if(r < 1.0){ 
        *lambdaT = Lambda/alpha;
        *lambdaI = r*(*lambdaT);
    }
    else { 
        *lambdaI = Lambda/alpha;
        *lambdaT = (*lambdaI)/r;
            //  printf("in get_lambdaIt...Lambda, r, lambdaI,T: %g %g %g %g \n", Lambda, r, *lambdaI, *lambdaT);
    }
} // end of get_LambdaIT_from_Lambda_r

void
get_lambdaIT_from_Lambda_r(double L_g, double Lambda, double r, double* lambdaI, double* lambdaT)
{
        // inverse of get_Lambda_alt
    double alpha = LAMBDA_DEF_FACTOR*L_g*L_g;
    if(r < 1.0){ 
        *lambdaT = Lambda/alpha;
        *lambdaI = r*(*lambdaT);
    }
    else { 
        *lambdaI = Lambda/alpha;
        *lambdaT = (*lambdaI)/r;
            //  printf("in get_lambdaIt...Lambda, r, lambdaI,T: %g %g %g %g \n", Lambda, r, *lambdaI, *lambdaT);
    }
} // end of get_LambdaIT_from_Lambda_r_dist

void
get_lambdaIT_from_Lambda_xi_old(int n_edges, double Lambda, double xi, double* lambdaI, double* lambdaT)
{
        // inverse of get_Lambda
    double alpha = LAMBDA_DEF_FACTOR*(double)NODISTLL;
    if(xi <= 1.0){ 
        *lambdaT = Lambda/alpha;
        *lambdaI = xi*(*lambdaT);
    }
    else { 
        *lambdaI = Lambda/alpha;
        *lambdaT = (*lambdaI)*(2.0 - xi);
            //  printf("in get_lambdaIt...Lambda, r, lambdaI,T: %g %g %g %g \n", Lambda, r, *lambdaI, *lambdaT);
    }
} // end of get_LambdaIT_from_Lambda_xi

void
get_lambdaIT_from_Lambda_xi(double L_g, double Lambda, double xi, double* lambdaI, double* lambdaT)
{
        // inverse of get_Lambda
    double alpha = LAMBDA_DEF_FACTOR*L_g*L_g;
    if(xi <= 1.0){ 
        *lambdaT = Lambda/alpha;
        *lambdaI = xi*(*lambdaT);
    }
    else { 
        *lambdaI = Lambda/alpha;
        *lambdaT = (*lambdaI)*(2.0 - xi);
            //  printf("in get_lambdaIt...Lambda, r, lambdaI,T: %g %g %g %g \n", Lambda, r, *lambdaI, *lambdaT);
    }
} // end of get_LambdaIT_from_Lambda_xi_dist 


Four_doubles update_lambdas(State* state, const Run_info_in* r_in, double exponent)
{
        // update lambdas
        // Mode 0; r = 2, no uniformization, Lambda fixed
        // Mode 1: r = 2, no uniformization, lambda simulated
        // Mode 2: Use uniformization, r = 2, lambdaI, lambdaT fixed
        // Mode 3: Use uniformization, lambdaI, lambdaT independent, simulated, 
        // Mode 4: Use uniformization, r = 2, Lambda fixed
        // Mode 5: Use uniformizatio, r = 2, Lambda simulated
        // Mode 6: Use uniformizatio, r, Lambda simulated
        // 0, 2, 4 should agree
        // 1, 5 should agree, also 3 if you only look at states near r = 2
        // 3, 6 should agree

  // exponent is 1 for cold chain, 1/T for hot chains.
    
        //   double lambdaI_prop, lambdaT_prop, Lambda_prop;
        //  double q_ratio, pi_ratio, p_accept;
        //  Lambdas* L = state->L;
        //   Lambdas L_prop = *(state->L); // initialize proposed lambdas to old lambdas

  // fprintf(stderr, "top of update_lambdas. state->L->r : %g\n", state->L->r);
  //  printf("in update_lambdas; exponent: %g\n", exponent);
  
    int NpM = (state->path->start->perm->n_edge);
    double L_g = state->path->start->perm->L_g;
    Four_doubles accepts = {0.0, 0.0, 0.0, 0.0};
    int iii, Nlambdaupdates = 4, Mode = r_in->Update_lambdas_mode, Nrupdates = 2;
    double Iavg, Tavg, Ai_avg, At_avg;
        //   double pi_ratio_1; // calculated using pi_ratio_gen, as a check
    int Nlambdaupdates_of_type;
        
    if(Mode == 3 || Mode == 7){
        Iavg = get_n_inv_avg(state->path); // average over steps on path of n_inv; depends on path, not on lambdas 
        Tavg = get_n_trans_avg(state->path); // average over steps on path of n_trans
        Ai_avg = get_A_i_avg(state->path);
        At_avg = get_A_t_avg(state->path);
    }
    
    Nlambdaupdates = 2*((Nlambdaupdates+1)/2); // if odd, go up to next even
    Nlambdaupdates_of_type = Nlambdaupdates;
    for(iii=0; iii<Nlambdaupdates; iii++){    
      if(r_in->Update_lambdas_mode == 7){ // in mode 7, update cold-chain lambdas alternating between mode 3 and 6.
            if(exponent == 1.0){ Mode = (iii%2 == 0)? 3: 6; Nlambdaupdates_of_type = Nlambdaupdates/2; }
            else { Mode = 3; }
        }
        
            // Modes 0, 1; no uniformization (this calculation assumes r = 2)
        if(Mode == 0){ return accepts;}// do nothing - leave lambdas as is
        else if(Mode == 1){ // no uniformization, just update Lambda with gibbs step
                // in this mode Lambda_max = lambda_max*NpMchoose2
                // state->L->Lambda = gamma_dev(state->path->length + 1, r_in->Lambda_max);
                // length + 2 because want prior for Lambda propto Lambda for comparison with indep lamdaI, lambdaT
            state->L->Lambda = gamma_dev(state->path->length + 2, r_in->Lambda_max);
            state->L->lambdaI = state->L->Lambda/(double)(NpM*(NpM-1)/2);
            state->L->lambdaT = state->L->lambdaI*0.5; 
        }
            // Modes 2, 3; use uniformization, lambdaI, lambdaT parametrization
        else if(Mode == 2){ return accepts;} // do nothing - leave lambdas as is
        else if(Mode == 3){
            accepts.a += update_lambdaI(state->path, state->L, LAMBDA_I_STEP_WIDTH_COEFF*sqrt((double)(state->path->length_i + 1))/Ai_avg,
                                        r_in->lambda_max, exponent)/(double)Nlambdaupdates_of_type;
            accepts.b += update_lambdaT(state->path, state->L, LAMBDA_T_STEP_WIDTH_COEFF*sqrt((double)(state->path->length_t + 1))/At_avg,
                                        r_in->lambda_max, exponent)/(double)Nlambdaupdates_of_type;     
        }
            // Modes 4, 5, 6; use uniformization, Lambda, r parametrization
        else if(Mode == 4){ return accepts; } // Lambda, r fixed
        else if(Mode >= 5){ // simulate Lambda
                // gibbs step for Lambda
                // arg is length+2 here because P(Lambda|X,r) ~ P(X|Lambda,r)P(Lambda) ~ e^(-Lambda)*Lambda^L * Lambda
                // = e^(-Lambda)*Lambda^(L+1).   The extra factor of Lambda comes from P(Lambda)

         /*    printf("length, length_a, lambdaI, lambdaT, : %i %i  %g %g \n", state->path->length, state->path->length_a, */
/*                   state->L->lambdaI, state->L->lambdaT); */
            state->L->Lambda = gamma_dev(state->path->length_a + 2, r_in->Lambda_max);
                //  printf("after gibbs step for Lambda\n");
            accepts.c += 1.0/(double)Nlambdaupdates_of_type;
            get_lambdaIT_from_Lambda_xi(L_g, state->L->Lambda, state->L->xi, &(state->L->lambdaI), &(state->L->lambdaT));
                // now *L is up to date
            if(Mode == 6){ // simulate xi also, rUPDATE_MODE is not 1
                int ii;

                for(ii=0; ii<Nrupdates; ii++){ // in params.h #define rUPDATE_MODE 1   , change to anything else to use xi parametrization
                    if(rUPDATE_MODE == 1) accepts.d += update_r(state->path, state->L)/(double)(Nlambdaupdates_of_type*Nrupdates);
                    else accepts.d += update_xi(state->path, state->L)/(double)(Nlambdaupdates_of_type*Nrupdates);         //  }
                } // loop over Nrupdates MH r updates
            } // end of if Mode == 6
        } // end of if Mode >= 5
    } // loop over lambda updates per call of update_lambda
        // getchar();
    assert(state->L->lambdaI >= 0.0 && state->L->lambdaI <= r_in->lambda_max);
    assert(state->L->lambdaT >= 0.0 && state->L->lambdaT <= r_in->lambda_max);
      
        //   printf("accepts.a,b,c,d: %8g %8g %8g %8g   %10g\n", accepts.a, accepts.b, accepts.c, accepts.d, exponent);
    return accepts;
} // end update_lambdas

void
print_Lambdas(Lambdas* L)
{
    printf("Lambda, r, xi, lambdaI, lambdaT: %g %g %g %g %g \n", L->Lambda, L->r, L->xi, L->lambdaI, L->lambdaT);
}


double
test_trans_sym(int size, long trans[][100])
{
    int i, j;
    double result = 0.0; long dof = (long)((size-1)*(size-2)/2);
    double x,y;
    for(i=1; i<size; i++){
        for(j=i+1; j<size; j++){
            x = (double)trans[i][j]; y = (double)trans[j][i];
            if(x+y > 0.5){
                result += (x-y)*(x-y)/((x+y)*(x+y));
            }
        }
    }
    return result/(double)dof;
}


State**
get_init_states(Permutation* p1, Permutation* p2, const Run_info_in* r_in, int n_short, int n_long, int max_path_length)
{
  State** the_state = (State**)chcalloc(r_in->N_chain, sizeof(State*)); // alloc array of pointers to the N_chain states
    int ii;
    double Epsilon_init;
    int flips[NMARKERMAX];
        //  double flip_probs[NMARKERMAX];
    int sum1=0, sum2=0, sum3 = 0, n_intermed = r_in->N_chain - (n_short + n_long), i_intermed;
        //   int length_t, length_f, length_a;
    static unsigned long firstseed;
    static int first = TRUE;
    unsigned long saveseed = rng_long;
    double sf_prob_long = r_in->Init_long_sf_param/sqrt((double)p1->n_mark);
    int Lest, Lestmin = MAXPATHLENGTH;
    Cycle_decomposition the_cd;
    double log_path_prob;
    Path* a_path;

    check_for_null_pointer((void*)the_state);
    
    get_CD(p1, p2, &the_cd); 
  
    if(first){
        firstseed = rng_long;
        first = FALSE;
    }
    else{ // If not updating paths, use same seed each time through, so if datasets are the same, paths will be also
            // seed the rng with firstseed  
    }
    if(!UPDATE_PATHS) seedrand(firstseed);
    fprintf(stdout, "Init path lengths: "); // getchar();

    
    for (ii = 0; ii<r_in->N_chain; ii++){ // gen initial paths, make n_short of them short, n_long long, the rest medium
            // initialize Imax etc.
      the_state[ii] = (State*)chcalloc(1, sizeof(State)); // allocate memory for State structure
      the_state[ii]->L = (Lambdas*)chcalloc(1, sizeof(Lambdas)); // allocate memory for Lambdas structure
        the_state[ii]->Imax = get_n_inv_max(p1);
        the_state[ii]->Imin = get_n_inv_min(p1);
        the_state[ii]->Tmax = get_n_trans_max(p1);
        the_state[ii]->Tmin = get_n_trans_min(p1);
        the_state[ii]->Mode = r_in->Update_lambdas_mode;
    }

    ii=0;
    while(ii<r_in->N_chain){
        
        i_intermed = ii - (n_short + n_long);
        if(r_in->Choose_signs > 0) { // unsigned case; good signs for short init paths, randomize signs for long init paths
           
            if (ii == 0) printf("Optimizing signs of p1 for short initial paths. \n");
            Lest = get_signs_short(p1, p2, r_in); // Lest could be 0
            printf("ii, Estimated length: %i %i \n", ii, Lest);
            Lestmin = (Lest < Lestmin)? Lest: Lestmin;
            
            if(ii < n_short){
                    // leave the signs "short"
            }
            else if(ii < n_short + n_long){
                (void)randomize_signs(p1->n_mark, sf_prob_long, flips); // flip some signs to increase distance
                flip_signs_by_position(p1, flips);
            }
            else { // intermediate length initial paths
                (void)randomize_signs(p1->n_mark, sf_prob_long*0.5, flips); // flip half as many signs to lengthen
                flip_signs_by_position(p1, flips);
            }
        }
        else { // signed case
            Lest = the_cd.n_mark - the_cd.n_chrom - the_cd.n_int_cycles + the_cd.c_g; 
            Lestmin = (Lest < Lestmin)? Lest: Lestmin;
        }
            // set Epsilon_init for short, medium, or long initial paths
        if(ii < n_short) {
            Epsilon_init = 0.5/(Lest+1); // so on average only 0.5 extra steps taken
        }
        else if(ii < n_short + n_long){
            Epsilon_init = (3.0/7.0)/(1.0 + sqrt((double)Lestmin)/r_in->Init_long_param);
        }
        else { // intermediate length initial paths
            Epsilon_init = 0.5*(0.5*Lest + (3.0/7.0)/(1.0 + sqrt((double)Lestmin)/r_in->Init_long_param)); // just avg of short and long cases for now
        }
        if(UPDATE_PATHS || ii == 0 /* || r_in->N_rev_fake > 0 */){

                // use INITLONGPATH_IPF*2*n_chrom as factor for preferring inversions over translocations in dd-0 or 1 steps in
                // initial path generation. This mean extra steps in long paths due to dd=0 or 1 (from larger epsilon)
                // are INITLONGPATH_IPF times as likely to be inverions as translocations. Extra length from sign flips is extra inversions
            if(TRUE || ii %2 == 0) {
                    //  print_perm(stdout, p1); print_perm(stdout, p2); //  getchar();
                printf("In get_init_states. before gen_path\n");
                a_path = gen_path(p1, p2, Epsilon_init, INITLONGPATH_IPF*2.0*(double)p1->n_chrom, max_path_length, &log_path_prob); // generate an initial main path
                printf("In get_init_states. after gen_path\n"); getchar();
            }
            else {
                a_path = gen_path(p2, p1, Epsilon_init, INITLONGPATH_IPF*2.0*(double)p1->n_chrom, max_path_length, &log_path_prob); // generate an initial main path
            }
            fflush(stdout);
            printf("a_path->length, max_path_length: %i  %i \n", a_path->length, max_path_length); fflush(stdout);
            if(a_path->length < max_path_length){
                the_state[ii]->path = a_path;
                printf("xxxxxxx ii, the_state[ii]->path->length: %i %i \n", ii, the_state[ii]->path->length );  fflush(stdout);
                printf("path length, length_t, length_i, L, L_t, L_i: %i %i %i   %i %i %i \n",
                       the_state[ii]->path->length, the_state[ii]->path->length_t,  the_state[ii]->path->length_i,
                       the_state[ii]->path->L, the_state[ii]->path->L_t, the_state[ii]->path->L_i);     
                if(ii < n_short /*  && ii%2 == 0 */){ // shorten /* half of */ the "short" init paths by removing loops
                        // shorten_path(the_state[ii]->path);
                }
                get_LLiLt(the_state[ii]->path);
                {
                    Cycle_decomposition the_cd;
                    get_CD(the_state[ii]->path->start->perm, the_state[ii]->path->end->perm, &the_cd); 
                    printf("distance: %i \n", the_cd.d);
                 printf("shortened path length, length_t, length_i, L, L_t, L_i: %i %i %i   %i %i %i \n",
                          the_state[ii]->path->length, the_state[ii]->path->length_t,  the_state[ii]->path->length_i,
                       the_state[ii]->path->L, the_state[ii]->path->L_t, the_state[ii]->path->L_i);     
                 printf("shortened path length: %i \n", the_state[ii]->path->length);  fflush(stdout); //getchar();
                }
                initialize_lambdas(the_state[ii], r_in);
                printf("ii: %i . the_state[ii]->L->lambdaI/T %g %g \n",
                       ii, the_state[ii]->L->lambdaI, the_state[ii]->L->lambdaT); 
                add_dummies(the_state[ii]->path, the_state[ii]->L);
                printf("L_a: %i \n", the_state[ii]->path->length_a);
                    // getchar();
                get_LLiLt(the_state[ii]->path);
                the_state[ii]->path->fission_factor = get_fission_factor(the_state[ii]->path);
                the_state[ii]->path->swap_end_factor = get_swap_end_factor(the_state[ii]->path);
                printf("in get_init_states, fission_factor: %g \n", the_state[ii]->path->fission_factor);  fflush(stdout);
                {
                    int check_path_lens, check_path_consist;
                    printf("lengths before, after shortening: %i ", the_state[ii]->path->length); fflush(stdout);
                    check_path_lens = (check_path_lengths(the_state[ii]->path) == TRUE);
                    check_path_consist = (check_path_consistency_multi(the_state[ii]->path) == 0); // checks that perms and reversals along path are consistent
                    printf("in get_init_state, after check_path_consistency...\n"); fflush(stdout);
                }
                ii++;   
            }
            else{ printf("in get_init_states, path too long, trying again. \n");  fflush(stdout); freepath(a_path); }// try again. generate another path - free this path
        }
        else{  // if not updating paths, just copy path to get init paths for other chains
                // so all chains have same path
            the_state[ii]->path = copy_path(the_state[0]->path);
            initialize_lambdas(the_state[ii], r_in);
        }
    }

    for(ii=0; ii<r_in->N_chain; ii++){
        
        if(ii < n_short) {sum1 += the_state[ii]->path->length;}
        else if(ii < n_short + n_long) {sum2 += the_state[ii]->path->length;}
        else {sum3 += the_state[ii]->path->length;}    
    
        assert(check_path_lengths(the_state[ii]->path) == TRUE);
    } //fprintf(stdout, "\n"); // end loop over chains     
  
    seedrand(saveseed); // restore rng

  /*   { */
/*         int ijk; */
/*         printf("Ntrans a,b: %g %g \n", (double)trans_a/r_in->N_chain, (double)trans_b/r_in->N_chain); */
/*          for(ijk=0; ijk<20; ijk++){ */
/*             fprintf(stdout, "%g %i %i \n", 0.05*(ijk+0.5), zeta_hist_a[ijk], zeta_hist_b[ijk]); */
/*          }fprintf(stdout, "\n"); getchar(); */
/*     } */
    return the_state;
} // end of get_init_states


Chain_set*
construct_chain_set(const Permutation* p1, const Permutation* p2, const Run_info_in* r_in)
{
        // construct an initialized set of N_chain * N_temperatures Chains
  Chain_set* the_set = (Chain_set*)chcalloc(1, sizeof(Chain_set)); // alloc Chain_set
    int i, j;
    double Temperature;

    printf("top of construct_chain_set\n");
    the_set->N_temperatures = r_in->N_temperatures;
    the_set->N_chain = r_in->N_chain;
        // chain is array of N_chain pointers to arrays of N_temperature pointers to Chains
    the_set->chain = (Chain***)chcalloc(r_in->N_chain, sizeof(Chain**));
        // now allocate the N_chain arrays of N_temperature arrays of Chains
    for(i=0; i<r_in->N_chain; i++){
      the_set->chain[i] = (Chain**)chcalloc(the_set->N_temperatures, sizeof(Chain*));
        Temperature = r_in->T_cold;
       printf("i,j,Temperature: %i %i %g \n", i, 0, Temperature); 
        the_set->chain[i][0] = construct_chain(p1, p2, r_in, i<NSHORT, Temperature);
        for(j=1; j<the_set->N_temperatures; j++){
             Temperature += (r_in->T_hot - r_in->T_cold)/(double)(r_in->N_temperatures - 1);
            printf("i,j,Temperature: %i %i %g \n", i, j, Temperature);
            the_set->chain[i][j] = construct_chain(p1, p2, r_in, i<NSHORT, Temperature);
        }
    }
        // allocate histograms
    printf("bottom of construct_chain_set\n");
    return the_set; 
} // end of construct_chain_set
        

Chain*
construct_chain(const Permutation* p1in, const Permutation* p2in, const Run_info_in* r_in, int make_short_path, double Temperature)
{
        // construct an initialized Chain
        //
  Chain* the_chain = (Chain*)chcalloc(1, sizeof(Chain)); // alloc array of pointers to the N_chain states
    State* the_state;
    Permutation* p1 = copy_perm(p1in);
    Permutation* p2 = copy_perm(p2in); 
    int i;
    double Epsilon_init;
    int flips[NMARKERMAX];
    double sf_prob_long = r_in->Init_long_sf_param/sqrt((double)p1->n_mark);
    int Lest, Lestmin = MAXPATHLENGTH;
    Cycle_decomposition the_cd;
    double log_path_prob;
    Path* a_path;
    int done = FALSE;

        //   printf("top of construct_chain\n");
    the_chain->temperature = Temperature;
//    printf("chain temperature: %g \n", the_chain->temperature); getchar();
    the_chain->N_steps_so_far = 0;
    the_chain->sumL = the_chain->sumLsqr = 0.0;
    the_chain->sumLi = the_chain->sumLisqr = 0.0;
    the_chain->sumLt = the_chain->sumLtsqr = 0.0;
    the_chain->sumlambdaI = the_chain->sumlambdaIsqr = 0.0;
    the_chain->sumlambdaT = the_chain->sumlambdaTsqr = 0.0;
    the_chain->L_min = MAXPATHLENGTH;
    the_chain->L_i_min = MAXPATHLENGTH;
    the_chain->L_t_min = MAXPATHLENGTH;
    the_chain->n_acc_path = 0.0;
    the_chain->n_acc_lambdaI = 0.0;
    the_chain->n_acc_lambdaT = 0.0;
    the_chain->n_acc_r = 0.0;
    the_chain->n_acc_Lambda = 0.0; 
    
    check_for_null_pointer((void*)the_chain);
    
    get_CD(p1, p2, &the_cd);
    
    the_state = (State*)chcalloc(1, sizeof(State)); // allocate memory for State structure
    the_state->L = (Lambdas*)chcalloc(1, sizeof(Lambdas)); // allocate memory for Lambdas structure
    the_state->Imax = get_n_inv_max(p1);
    the_state->Imin = get_n_inv_min(p1);
    the_state->Tmax = get_n_trans_max(p1);
    the_state->Tmin = get_n_trans_min(p1);
    the_state->Mode = r_in->Update_lambdas_mode;
    
    for(i=0; !done; i++){ // loop until acceptable (i.e. short enough) path is found 
            // fiddle with signs (in unsigned case) for short or long paths
        if(r_in->Choose_signs > 0) { // unsigned case; good signs for short init paths, randomize signs for long init paths
           
            printf("Optimizing signs of p1 for short initial paths. \n"); 
            Lest = get_signs_short(p1, p2, r_in); // Lest could be 0
            printf("Estimated length: %i \n", Lest);
            
            if (!make_short_path) {
                (void)randomize_signs(p1->n_mark, sf_prob_long, flips); // flip some signs to increase distance
                flip_signs_by_position(p1, flips);
            }
        }
        else { // signed case
            Lest = the_cd.n_mark - the_cd.n_chrom - the_cd.n_int_cycles + the_cd.c_g;
        }
        Lestmin = (Lest < Lestmin)? Lest: Lestmin;
         
            // set Epsilon_init for short, or long initial paths
        if(make_short_path) {
            Epsilon_init = 0.5/(Lest+1); // so on average only 0.5 extra steps taken
        }
        else {
            Epsilon_init = (3.0/7.0)/(1.0 + sqrt((double)Lestmin)/r_in->Init_long_param);
        }

        printf("in construct_chain. after Epsilon_init = ...; Epsilon_init: %g \n", Epsilon_init);
            // use INITLONGPATH_IPF*2*n_chrom as factor for preferring inversions over translocations in dd-0 or 1 steps in
            // initial path generation. This mean extra steps in long paths due to dd=0 or 1 (from larger epsilon)
            // are INITLONGPATH_IPF times as likely to be inverions as translocations. Extra length from sign flips is extra inversions
                   
        a_path = gen_path(p1, p2, Epsilon_init, INITLONGPATH_IPF*2.0*(double)p1->n_chrom, r_in->MC3_max_pathlength, &log_path_prob); // generate an initial main path
        printf("in construct_chain. after gen_path = ...; length: %i %i \n", a_path->length, make_short_path);
        printf("check path revs ordered. %i \n", check_path_revs_ordered(a_path));
        if(make_short_path){
                  a_path = shorten_path(a_path); // shorten the "short" init paths by removing loops
            printf("in construct_chain. after shorten_path = ...; length: %i \n", a_path->length);
        }
        if(a_path->length < r_in->MC3_max_pathlength){
            the_state->path = a_path;
	    fprintf(stderr, "L_a: %ld\n", a_path->length); //getchar();
               
            get_LLiLt(the_state->path); // get L, L_i, L_t, i.e. lengths not counting within equivalence class transitions
            {
                Cycle_decomposition the_cd;
                get_CD(the_state->path->start->perm, the_state->path->end->perm, &the_cd); 
                printf("distance: %i \n", the_cd.d);
                printf("shortened path length, length_t, length_i, L, L_t, L_i: %i %i %i   %i %i %i \n",
                       the_state->path->length, the_state->path->length_t,  the_state->path->length_i,
                       the_state->path->L, the_state->path->L_t, the_state->path->L_i);     
                printf("shortened path length: %i \n", the_state->path->length);  fflush(stdout); //getchar();
            }
            initialize_lambdas(the_state, r_in);
          
            add_dummies(the_state->path, the_state->L);
            printf("L_a: %i \n", the_state->path->length_a);  //getchar();
                   
            get_LLiLt(the_state->path);
            the_state->path->fission_factor = get_fission_factor(the_state->path);
            the_state->path->swap_end_factor = get_swap_end_factor(the_state->path);
            printf("in construct_chain, fission_factor: %g \n", the_state->path->fission_factor);  fflush(stdout);
            {
                int check_path_lens, check_path_consist;
                printf("lengths before, after shortening: %i ", the_state->path->length); fflush(stdout);
                check_path_lens = (check_path_lengths(the_state->path) == TRUE);
                check_path_consist = (check_path_consistency_multi(the_state->path) == 0); // checks that perms and reversals along path are consistent
                printf("in get_init_state, after check_path_consistency...\n"); fflush(stdout);
            }
            done = TRUE;
        }
        else{
            printf("in construct_chain, path too long, trying again. \n"); freepath(a_path);
            Epsilon_init *= 0.95; sf_prob_long *= 0.95; // so next path should be a bit shorter
            if(i >= 20){
                printf("In construct_chain. generated >20 paths and all too long. Exiting\n"); process_error();
            }
        }// try again. generate another path - free this path

    }
    the_chain->state = the_state;

    // fprintf(stderr, "In construct_chain: the_state->L->r: %g\n", the_state->L->r); getchar();

        // allocate histograms
     printf("bottom of construct_chain\n");
    return the_chain;
} // end of construct_chain

void
update_chain(Chain* the_chain, const Run_info_in* r_in)
{
    State* the_state = the_chain->state;
    Path* the_path = the_state->path;
    Four_doubles lambda_accepts;
    
    if(FIXEDLAMBDAr){ 
        double lamI, lamT;
        double L_g = the_state->path->start->perm->L_g;
                    
        the_state->L->Lambda = FIXED_LAMBDA;
        the_state->L->r = FIXED_r;
        get_lambdaIT_from_Lambda_r(L_g, the_state->L->Lambda, the_state->L->r, &lamI, &lamT);
        the_state->L->lambdaI = lamI; the_state->L->lambdaT = lamT;
        the_state->L->xi = get_xi(the_state->L->lambdaI, the_state->L->lambdaT);
    }
    else{
        lambda_accepts = update_lambdas(the_state, r_in, 1.0/the_chain->temperature);
    }
               
    the_chain->n_acc_lambdaI += lambda_accepts.a;
    the_chain->n_acc_lambdaT += lambda_accepts.b;
    the_chain->n_acc_r += lambda_accepts.d;
    the_chain->n_acc_Lambda += lambda_accepts.c;
                
    if(UPDATE_PATHS){
        if(the_path->length > r_in->MC3_max_pathlength){ printf("MC3_max_pathlength, the_path->length: %i %i  \n", r_in->MC3_max_pathlength, the_path->length);}
        the_chain->n_acc_path += update_path(the_state, r_in, 1.0/the_chain->temperature, UPDATE_TO_END); //, n_props, n_accepts, totalLprop_hist, totalLacc_hist);
            //     printf("after update_path. ii, iii, accepted: %i  %i  %i \n", ii, iii, accepted); 
            //   update_path_order(the_state->path);
    }
    else{
        the_chain->n_acc_path += 1; 
    }
    the_path = the_state->path;
    the_chain->N_steps_so_far++;
    the_chain->sumL += the_path->length; the_chain->sumLsqr += the_path->length*the_path->length;
    the_chain->sumLi += the_path->length_i; the_chain->sumLisqr += the_path->length_i*the_path->length_i;
    the_chain->sumLt += the_path->length_t; the_chain->sumLtsqr += the_path->length_t*the_path->length_t;
    the_chain->sumlambdaI += the_state->L->lambdaI; the_chain->sumlambdaIsqr += the_state->L->lambdaI*the_state->L->lambdaI;
    the_chain->sumlambdaT += the_state->L->lambdaT; the_chain->sumlambdaTsqr += the_state->L->lambdaT*the_state->L->lambdaT;

    the_chain->W_L = within_chain_variance(the_chain->N_steps_so_far, the_chain->sumL, the_chain->sumLsqr); // get the within-chain variance
    the_chain->W_Li = within_chain_variance(the_chain->N_steps_so_far, the_chain->sumLi, the_chain->sumLisqr);
    the_chain->W_Lt = within_chain_variance(the_chain->N_steps_so_far, the_chain->sumLt, the_chain->sumLtsqr);
    the_chain->W_lambdaI = within_chain_variance(the_chain->N_steps_so_far, the_chain->sumlambdaI, the_chain->sumlambdaIsqr);
    the_chain->W_lambdaT = within_chain_variance(the_chain->N_steps_so_far, the_chain->sumlambdaT, the_chain->sumlambdaTsqr);

    if(the_path->L < the_chain->L_min) the_chain->L_min = the_path->L;
    if(the_path->L_i < the_chain->L_i_min) the_chain->L_i_min = the_path->L_i;
   if(the_path->L_t < the_chain->L_t_min) the_chain->L_t_min = the_path->L_t; 
  
} // end of update_chain

inline double
within_chain_variance(int N, double sumX, double sumXsqr)
{
    double avgX = sumX/(double)N;
    double avgXsqr = sumXsqr/(double)N;
    return avgXsqr - avgX*avgX;
}

void
chain_set_abw(Chain_set* chain_set)
{
        // update the values of Avg, B and W (average, and between and within chain variances) for L, Li, etc.
        // for each temperature by looking at the chains

    int i, j, N_ch = chain_set->N_chain, N_T = chain_set->N_temperatures;
    Chain* the_chain;
    double WL, AveL, AveAveL = 0.0, sumAveLsqr = 0.0;
    double WLi, AveLi, AveAveLi = 0.0, sumAveLisqr = 0.0;
    double WLt, AveLt, AveAveLt = 0.0, sumAveLtsqr = 0.0;
    double WlambdaI, AvelambdaI, AveAvelambdaI = 0.0, sumAvelambdaIsqr = 0.0;
    double WlambdaT, AvelambdaT, AveAvelambdaT = 0.0, sumAvelambdaTsqr = 0.0;  
    
    for(i=0; i<N_T; i++){
        WL = WLi = WLt = WlambdaI = WlambdaT = 0.0;
        AveAveL =  AveAveLi = AveAveLt = AveAvelambdaI = AveAvelambdaT = 0.0;
        sumAveLsqr = sumAveLisqr = sumAveLtsqr = sumAvelambdaIsqr = sumAvelambdaTsqr = 0.0;
        for(j=0; j<N_ch; j++){
            the_chain = chain_set->chain[j][i];
            WL += the_chain->W_L;
            AveL = the_chain->sumL/(double)the_chain->N_steps_so_far;
            AveAveL += AveL; sumAveLsqr += AveL*AveL;

            WLi += the_chain->W_Li;
            AveLi = the_chain->sumLi/(double)the_chain->N_steps_so_far;
            AveAveLi += AveLi; sumAveLisqr += AveLi*AveLi;

            WLt += the_chain->W_Lt;
            AveLt = the_chain->sumLt/(double)the_chain->N_steps_so_far;
            AveAveLt += AveLt; sumAveLtsqr += AveLt*AveLt;

            WlambdaI += the_chain->W_lambdaI;
            AvelambdaI = the_chain->sumlambdaI/(double)the_chain->N_steps_so_far;
            AveAvelambdaI += AvelambdaI; sumAvelambdaIsqr += AvelambdaI*AvelambdaI;

            WlambdaT += the_chain->W_lambdaT;
            AvelambdaT = the_chain->sumlambdaT/(double)the_chain->N_steps_so_far;
                // printf("AvelambdaT: %g \n", AvelambdaT);
            AveAvelambdaT += AvelambdaT; sumAvelambdaTsqr += AvelambdaT*AvelambdaT;            
        }
        
        chain_set->W_L[i] = WL/(double)N_ch; 
        AveAveL /= (double)N_ch; chain_set->Avg_L[i] = AveAveL;
        chain_set->B_L[i] = sumAveLsqr/(double)N_ch - AveAveL*AveAveL;

        chain_set->W_Li[i] = WLi/(double)N_ch; 
        AveAveLi /= (double)N_ch; chain_set->Avg_Li[i] = AveAveLi;
        chain_set->B_Li[i] = sumAveLisqr/(double)N_ch - AveAveLi*AveAveLi;

        chain_set->W_Lt[i] = WLt/(double)N_ch; 
        AveAveLt /= (double)N_ch; chain_set->Avg_Lt[i] = AveAveLt;
        chain_set->B_Lt[i] = sumAveLtsqr/(double)N_ch - AveAveLt*AveAveLt;

        chain_set->W_lambdaI[i] = WlambdaI/(double)N_ch; 
        AveAvelambdaI /= (double)N_ch; chain_set->Avg_lambdaI[i] = AveAvelambdaI;
        chain_set->B_lambdaI[i] = sumAvelambdaIsqr/(double)N_ch - AveAvelambdaI*AveAvelambdaI;

        chain_set->W_lambdaT[i] = WlambdaT/(double)N_ch; 
        AveAvelambdaT /= (double)N_ch; chain_set->Avg_lambdaT[i] = AveAvelambdaT;
        
        chain_set->B_lambdaT[i] = sumAvelambdaTsqr/(double)N_ch - AveAvelambdaT*AveAvelambdaT;
        
     /*    printf("i, B,W lambdaI, B,W lambdaT: %i  %g %g %g %g \n", i, chain_set->B_lambdaI[i], chain_set->W_lambdaI[i], */
/*                chain_set->B_lambdaT[i],  chain_set->W_lambdaT[i]); */
    }       
}

int 
X_converged(int index, const Run_info_in* r_in, double** cumeSteps, double** cumeX, double** cumeXsqr, struct XBW* the_XBW)
{  // returns true if set of chains has converged
 
  int ii;
  int n_half, n_whole, half_index, has_converged = FALSE;     
  double sum_X_2[NCHAINMAX] = {0.0}, sum_Xsqr_2[NCHAINMAX] = {0.0};
  double Xave[NCHAINMAX], S_sqr[NCHAINMAX]; 
  double Xave_2[NCHAINMAX], S_sqr_2[NCHAINMAX];
  double W=0.0, B=0.0, Xaverage=0.0, srh;
  double W_2=0.0, B_2=0.0, Xaverage_2=0.0, srh_2; 

  if(index < 2*N_DOUBLE){
      half_index = index/2;
      if(2*half_index == index || index <= 2) return FALSE;
  }
  else half_index = index - N_DOUBLE;
  
  for (ii=0; ii<r_in->N_chain; ii++){
      n_whole = (int)cumeSteps[index][ii];
      n_half = (int)cumeSteps[half_index][ii];
          // if(n_whole % 2 != 0) {printf("in X_converged, n_whole is not even. n_whole: %i \n", n_whole); return FALSE;} // 
      if (n_whole != 2*n_half) { 
          printf("in X_converged, n_whole != 2*n_half\n");
          process_error();
      } 
      sum_X_2[ii] = cumeX[index][ii] - cumeX[half_index][ii];
      sum_Xsqr_2[ii] = cumeXsqr[index][ii] - cumeXsqr[half_index][ii];

      Xave[ii] = cumeX[index][ii]/n_whole; 
      Xaverage += Xave[ii]/(double)r_in->N_chain; 
      S_sqr[ii] = cumeXsqr[index][ii]/n_whole - Xave[ii]*Xave[ii];
      if(S_sqr[ii] < -1.0e-7) {
          printf("S_sqr[ii]: %g, index: %i, ii: %i, cumeXsqr[index][ii]/(n_whole): %g, Xave[ii]: %g \n",
                 S_sqr[ii], index, ii, cumeXsqr[index][ii]/(n_whole), Xave[ii]);
          printf("index, ii, n_whole, cumeXsqr[index][ii], cumeX[index][ii]: %i %i  %i  %g %g \n", index, ii, n_whole, cumeXsqr[index][ii], cumeX[index][ii]);
          process_error();
      }

      B += Xave[ii]*Xave[ii]/(double)r_in->N_chain; 
      W += S_sqr[ii]/(double)r_in->N_chain; 

      Xave_2[ii] = sum_X_2[ii]/(n_half); 
      Xaverage_2 += Xave_2[ii]/(double)r_in->N_chain; 
      S_sqr_2[ii] = sum_Xsqr_2[ii]/(n_half) - Xave_2[ii]*Xave_2[ii]; 
      B_2 += Xave_2[ii]*Xave_2[ii]/(double)r_in->N_chain; 
      W_2 += S_sqr_2[ii]/(double)r_in->N_chain;
  }           
 
  B -= Xaverage*Xaverage; 
  B *= (double)r_in->N_chain/(double)(r_in->N_chain - 1); 
  B_2 -= Xaverage_2*Xaverage_2; 
  B_2 *= (double)r_in->N_chain/(double)(r_in->N_chain - 1); 
 
 
  if(!((W>0.0) && (W_2>0.0))){printf("In X_converged. Warning: W, W_2: %f %f \n", W, W_2);/*  process_error(); */} 
  srh = (W>0.0)? sqrt((B+W)/W): 10000.0; 
  srh_2 = (W_2>0.0)? sqrt((B_2+W_2)/W_2): 10000.0;

  the_XBW->X = Xaverage; the_XBW->B = B; the_XBW->W = W;
  the_XBW->X2 = Xaverage_2; the_XBW->B2 = B_2; the_XBW->W2 = W_2;
//  printf("Xaverage, B, W: %g %g %g \n", Xaverage, B, W);
  
  if ((srh < r_in->Conv_srh) && (srh_2 < r_in->Conv_srh_half2)){ 
      has_converged = TRUE; 
  }
  
  return has_converged; 
}  // end of X_converged
 
void check_run_params(const Run_info_in* r_in) 
{ 
  // checks parameter values against max values allowed
  if ((r_in->N_chain > NCHAINMAX) || (r_in->N_chain < NCHAINMIN)){ 
    printf("Warning: N_chain value out of range. \n"); 
    printf("Allowed range %i -> %i. N_chain: %i.\n Hit any key to continue \n", 
	   NCHAINMIN, NCHAINMAX, r_in->N_chain); 
    getchar();
  } 
  if ((r_in->N_data > NDATAMAX) || (r_in->N_data < NDATAMIN)){ 
    printf("Warning: N_data value out of range. \n"); 
    printf("Allowed range %i -> %i. N_data: %i.\n Hit any key to continue \n", 
	   NDATAMIN, NDATAMAX, r_in->N_data); 
    getchar();
  } 
  if ((r_in->Max_runlength > MAXRUNLENGTH) || (r_in->Max_runlength < MINRUNLENGTH)){ 
    printf("Warning: Max_runlength value out of range. \n"); 
    printf("Allowed range %i -> %i. Max_runlength: %i.\n Hit any key to continue \n", 
	   MINRUNLENGTH, MAXRUNLENGTH, r_in->Max_runlength);
    getchar(); 
  } 
} // end of check_run_params 



void input_run_params(FILE* fp, Run_info_in* r_in) 
{
    fscanf(fp, "%*s%s %*s%i %*s%i %*s%lf %*s%i %*s%i",  
           r_in->Genome_data_file,
           &r_in->Use_distances, 
           &r_in->N_rev_fake,
           &r_in->lambda_ratio_fake, 
           &r_in->N_chain,
           &r_in->N_data);
    fscanf(fp, "%*s%i  %*s%lf %*s%lf",
           &r_in->Update_lambdas_mode,
           &r_in->Lambda,
           &r_in->lambda_max);
    r_in->Lambda_max = UNKNOWN; // gets set in grit.c
    
    fscanf(fp, "%*s%lf  %*s%lf  %*s%lf ", 
               //    &r_in->Epsilon_init_long, &r_in->Epsilon_init_short,
           &r_in->Epsilon,
           &r_in->Init_long_param, &r_in->Init_long_sf_param);

   
     fscanf(fp, "%*s%i  %*s%lf  %*s%lf  %*s%i %*s%i ", 
               //    &r_in->Epsilon_init_long, &r_in->Epsilon_init_short,
           &r_in->N_temperatures, &r_in->T_cold, &r_in->T_hot, 
           &r_in->Do_MC3_swap, &r_in->MC3_max_pathlength); 

    fscanf(fp, "%*s%lf %*s%lf  %*s%lf %*s%lf  %*s%lf  %*s%lf %*s%lf   %*s%i %*s%i %*s%i %*s%i",
           &r_in->Conv_srh,
           &r_in->Conv_srh_half2,
           &r_in->Conv_Ldist,
           &r_in->Conv_Lambdadist,
           &r_in->Conv_rdist,
           &r_in->Conv_lambdaIdist,
           &r_in->Conv_lambdaTdist,
           &r_in->Burn_in,
           &r_in->Stop_at,
           &r_in->Min_burn_in_length,
           &r_in->Max_runlength);

    fscanf(fp, "%*s%i %*s%i  %*s%lf %*s%lf  %*s%i  %*s%lf",
           &r_in->Rand_seed,
           &r_in->L_prop_dist_type,
           &r_in->Alpha,
           &r_in->Xi,
           &r_in->Choose_signs,
               //          &r_in->Flip_some_sign_epsilon,
           &r_in->Flip_one_sign_epsilon); 
} // end of input_run_paramss


void output_run_params(FILE* fp, const Run_info_in* r_in) 
{
  char fc = '#';
    fprintf(fp, "%c %-24s%s\n", fc, "Genome_data_file", r_in->Genome_data_file);
    fprintf(fp, "%c %-24s%i\n %c %-24s%i\n %c %-24s%g\n %c %-24s%i\n %c %-24s%i\n#\n",
            fc, "Use_distances", r_in->Use_distances,
            fc, "N_rev_fake", r_in->N_rev_fake,
            fc, "lambda_ratio_fake", r_in->lambda_ratio_fake, 
            fc, "N_chain", r_in->N_chain,
            fc, "N_data", r_in->N_data);
    
    fprintf(fp, "%c %-24s%i\n%c %-24s%g\n%c %-24s%g\n%c %-24s%g\n#\n",
            fc, "Update_lambdas_mode", r_in->Update_lambdas_mode,
            fc, "Lambda", r_in->Lambda,
            fc, "lambda_max", r_in->lambda_max, 
            fc, "Lambda_max", r_in->Lambda_max);
    
    fprintf(fp, "%c %-24s%g\n %c %-24s%g\n %c %-24s%g\n#\n",
                //  "Epsilon_inits", r_in->Epsilon_init_long, r_in->Epsilon_init_short,
            fc, "Epsilon", r_in->Epsilon,
            fc, "Init_long_param", r_in->Init_long_param,
	    fc, "Init_long_sf_param", r_in->Init_long_sf_param);

    fprintf(fp, "%c %-24s%i\n%c %-24s%g\n%c %-24s%g\n%c %-24s%i\n%c %-24s%i\n#\n",
            fc, "N_temperatures",r_in->N_temperatures, fc, "T_cold", r_in->T_cold, fc, "T_hot", r_in->T_hot, 
            fc, "Do_MC3_swap", r_in->Do_MC3_swap, fc, "MC3_max_pathlength", r_in->MC3_max_pathlength); 

    fprintf(fp, "%c %-24s%g\n%c %-24s%g\n%c %-24s%g\n%c %-24s%g\n%c %-24s%g\n%c %-24s%g\n%c %-24s%g\n%c %-24s%i\n%c %-24s%i\n%c %-24s%i\n%c %-24s%i\n#\n",
            fc, "Conv_srh", r_in->Conv_srh,
            fc, "Conv_srh_half2", r_in->Conv_srh_half2,
            fc, "Conv_Ldist", r_in->Conv_Ldist,
            fc, "Conv_Lambdadist", r_in->Conv_Lambdadist,
            fc, "Conv_rdist", r_in->Conv_rdist,
            fc, "Conv_lambdaIdist", r_in->Conv_lambdaIdist,
            fc, "Conv_lambdaTdist", r_in->Conv_lambdaTdist,
            fc, "Burn_in", r_in->Burn_in,
            fc, "Stop_at", r_in->Stop_at,
            fc, "Min_burn_in_length ", r_in->Min_burn_in_length,
            fc, "Max_runlength", r_in->Max_runlength);

    fprintf(fp,  "%c %-24s%i\n%c %-24s%i\n%c %-24s%g\n%c %-24s%g\n%c%-24s%i\n%c %-24s%g\n#\n",
            fc, "Rand_seed", (int)r_in->Rand_seed,
                //"Path_seed_inc", (int)r_in->Path_seed_inc,
                //"MC_seed", (int)r_in->MC_seed,
            fc, "L_prop_dist_type", r_in->L_prop_dist_type,
            fc, "Alpha", r_in->Alpha,
            fc, "Xi", r_in->Xi,
            fc, "Choose_signs", r_in->Choose_signs,
                //          fc, "Flip_some_sign_epsilon", r_in->Flip_some_sign_epsilon,
            fc, "Flip_one_sign_epsilon", r_in->Flip_one_sign_epsilon);
#ifdef NDEBUG     // debugging is off
    fprintf(fp, "# ,NDEBUG defined. Debugging off. \n");
#else
    fprintf(fp, "# ,NDEBUG not defined. Debugging (i.e. assert) is on. \n");
#endif
    fprintf(fp, "# ,UPDATE_WHOLE_PATH: %i \n", UPDATE_WHOLE_PATH);
    fprintf(fp, "# ,USE2EPSILONS: %i \n", USE2EPSILONS); 
    fprintf(fp, "# ,EPSILON_BIG_FACTOR: %g \n", EPSILON_BIG_FACTOR); 
    fprintf(fp, "# ,BIGEPSILONFRACTION: %g \n", BIGEPSILONFRACTION);        
    fprintf(fp, "# ,UPDATE_PATHS: %i \n", UPDATE_PATHS);
    fprintf(fp, "# ,MAXPATHLENGTH: %i \n", MAXPATHLENGTH);
    fprintf(fp, "# ,MAXSHORTPATHLENGTH: %i \n", MAXSHORTPATHLENGTH);
    fprintf(fp, "# ,Using log(ratios).\n#\n");
 // *****************************************   
    fprintf(fp, "# ,LAMBDA_DEF_FACTOR: %g \n", LAMBDA_DEF_FACTOR);
    fprintf(fp, "# ,UPDATE_lambdaIT_MODE: %i \n", UPDATE_lambdaIT_MODE);
    fprintf(fp, "# ,LAMBDA_I_STEP_WIDTH_COEFF: %g \n", LAMBDA_I_STEP_WIDTH_COEFF);
    fprintf(fp, "# ,LAMBDA_T_STEP_WIDTH_COEFF: %g \n", LAMBDA_T_STEP_WIDTH_COEFF);
    fprintf(fp, "# ,rUPDATE_MODE: %i \n", rUPDATE_MODE);
    fprintf(fp, "# ,r_STEP_WIDTH: %g \n", r_STEP_WIDTH); 
    fprintf(fp, "# ,XI_WIDTH: %g \n#\n", XI_WIDTH);
    fprintf(fp, "# ,PREFERDD0SEDS: %i \n", PREFERDD0SEDS); 
    fprintf(fp, "# ,PDD0SEDS: %g \n", PDD0SEDS);
    fprintf(fp, "# ,PREFERDD0SEDS_SIGNFLIP: %i \n", PREFERDD0SEDS_SIGNFLIP);
    fprintf(fp, "# ,PDD0SEDS_SIGNFLIP: %g \n#\n", PDD0SEDS_SIGNFLIP);
    fprintf(fp, "# ,PREFER_INVERSIONS: %i \n",  PREFER_INVERSIONS);
        //   fprintf(fp, "# ,INVERSION_PREFER_FACTOR: %s \n", INVERSION_PREFER_FACTOR);
    fprintf(fp, "# ,USE_IPF_LAMBDA_RATIO: %i \n", USE_IPF_LAMBDA_RATIO);
    fprintf(fp, "# ,IPF_EXPONENT: %g \n", IPF_EXPONENT);
    fprintf(fp, "# ,ddm1_prefer_inversions: %i \n", ddm1_prefer_inversions);
    fprintf(fp, "# ,DDM1_IT_PROB_RATIO: %g \n", DDM1_IT_PROB_RATIO);
        //  fprintf(fp, "# ,DDM1_IT_PROB_RATIO: %i \n", DDM1_IT_PROB_RATIO);
    fprintf(fp, "# ,COUNT_WITHIN_EQUIVALENCE_CLASS_STEPS: %i \n", COUNT_WITHIN_EQUIVALENCE_CLASS_STEPS);
    fprintf(fp, "# ,LTHEAT: %i \n", LTHEAT);
    fprintf(fp, "# ,LTHEATFACTOR: %g \n", LTHEATFACTOR);
    fprintf(fp, "# ,LTHEATKNEE: %i \n", LTHEATKNEE);
    fprintf(fp, "# ,LTSHARP: %g \n", LTSHARP);
    fprintf(fp, "# ,INITLONGPATH_IPF: %g \n", INITLONGPATH_IPF);
    fprintf(fp, "\n");
}// end output_run_params


void
output2(int nsteps, struct XBW ABW, struct XBW LBW, struct XBW lambdaBW)
{
    fprintf(fp_out2, "%i %9g %9g %9g %9g %9g %9g \n%i %9g %9g %9g %9g %9g %9g \n%i %9g %9g %9g %9g %9g %9g \n", 
            nsteps, ABW.X, ABW.B, ABW.W, ABW.X2, ABW.B2, ABW.W2, 
            nsteps, LBW.X, LBW.B, LBW.W, LBW.X2, LBW.B2, LBW.W2,
            nsteps, lambdaBW.X, lambdaBW.B, lambdaBW.W, lambdaBW.X2, lambdaBW.B2, lambdaBW.W2);
}

void
output4(FILE* fp, Path_tree* path_tree, Histogram* L_hist, Histogram* Lambda_hist,
        Histogram* L_i_hist, Histogram* lambdaI_hist, Histogram* L_t_hist, Histogram* lambdaT_hist)
{
    print_path_tree(fp, path_tree, TRUE);
    print_histogram(fp, L_hist, 1, 1.0);
    print_histogram(fp, Lambda_hist, 3, 1.0); 
    print_histogram(fp, L_i_hist, 1, 1.0);
    print_histogram(fp, lambdaI_hist, 3, 1.0);
    print_histogram(fp, L_t_hist, 1, 1.0);
    print_histogram(fp, lambdaT_hist, 3, 1.0);
}

void output5_old(FILE* fp, int steps_so_far, int i, double** cumeL, int N_chain)
{        // output chain_averages stuff
    int ii, half_index;
    double L2;
     
    fprintf(fp, "%i ", steps_so_far);
    for(ii=0; ii<N_chain; ii++){
        half_index = i - N_DOUBLE;
        if(half_index >= N_DOUBLE){
            L2 = (cumeL[i][ii]-cumeL[half_index][ii])*(2.0/steps_so_far);
        }
        else if(i == 0){
            L2 = cumeL[i][ii];
        }
        else {
            half_index = (i+1)/2 - 1;         
            L2 = (cumeL[i][ii]-cumeL[half_index][ii])/(NAVG0*(i-half_index));
        }
        fprintf(fp, "%g %g ", cumeL[i][ii]/steps_so_far, L2);
    }fprintf(fp, "\n");
}

void output5(FILE* fp, int steps_so_far, int i, double** cumeL, double** cumeLambda, int* Lengths, int N_chain)
{        // output chain_averages stuff
    int ii, half_index;
    double L2;
     
    fprintf(fp, "%i  ", steps_so_far);
    for(ii=0; ii<N_chain; ii++){
        half_index = i - N_DOUBLE;
        if(half_index >= N_DOUBLE){
            L2 = (cumeL[i][ii]-cumeL[half_index][ii])*(2.0/steps_so_far);
        }
        else if(i == 0){
            L2 = cumeL[i][ii];
        }
        else {
            half_index = (i+1)/2 - 1;         
            L2 = (cumeL[i][ii]-cumeL[half_index][ii])/(NAVG0*(i-half_index));
        }
        fprintf(fp, "%i %g %g ", Lengths[ii], cumeL[i][ii]/steps_so_far, cumeLambda[i][ii]/steps_so_far);
    }fprintf(fp, "\n");
}

void
output7(FILE* fp, int Nchains, Histogram** L_hist, Histogram** lambda_hist)
{ 
        //output_run_params(fp);
        //print_histograms(fp, Nchains, L_hist, 1);
        //print_histograms(fp, Nchains, lambda_hist, 1);
  
  print_histograms(fp, Nchains, L_hist, 1, 1.0);
  print_histograms(fp, Nchains, lambda_hist, 3, 1.0);
 
  fflush(fp);
}


double
dist_converge(int N_chain, Histogram** hist)
{
    Stats* hist_stats;
    Histogram* the_hist = hist[0];
    int i, j, nbins = the_hist->nbins;
    double sum_means = 0.0, sum_stddevs = 0.0;

    hist_stats = (Stats*)chcalloc(nbins, sizeof(Stats));
    check_for_null_pointer((void*)hist_stats);

    for(j=0; j<nbins; j++){
        initialize_stats(hist_stats+j);
    }

        // get the mean (over chains) and stddev for each bin
    for(i=0; i< N_chain; i++){
        the_hist = hist[i];
        for(j=0; j<nbins; j++){
            insert_in_stats(hist_stats + j, the_hist->histogram[j], 1.0);
        }
    }

    for(j=0; j<nbins; j++){
        sum_means += hist_stats[j].mean;
        sum_stddevs += hist_stats[j].stddev;
    }

    free(hist_stats);
    return (sum_stddevs/sum_means);
}

double
check_path_freq_vs_pi(Path_tree* tree)
{

    Path_tree_node* node = tree->start;
    double sum_n = 0.0, sum_p = 0.0, result = 0.0;

    if(node == NULL) return 0.0;

        // get sums
    do{
        sum_n += node->count;
        sum_p += node->pi_rel;
        node = node->right;
    }while(node != NULL); //

    node = tree->start;
    do{
        result += (node->count/sum_n - (node->pi_rel)/sum_p)*
            (node->count/sum_n - node->pi_rel/sum_p)*(sum_p/node->pi_rel);
       node = node->right; 
    }while(node != NULL); //
    return result;
}  // end of check_path_freq_vs_pi

double
check_path_freq_vs_pi_old(Path_tree* tree)
{

    Path_tree_node* node = tree->start;
    double sum_n = 0.0, sum_p = 0.0, result = 0.0;

    if(node == NULL) return 0.0;

        // get sums
    do{
        sum_n += node->count;
        sum_p += node->pi_rel_old;
        node = node->right;
    }while(node != NULL); //

    node = tree->start;
    do{
        
        result += (node->count/sum_n - node->pi_rel_old/sum_p)*
            (node->count/sum_n - node->pi_rel_old/sum_p)*(sum_p/node->pi_rel_old);
       node = node->right; 
    }while(node != NULL); //
    return result;
}  // end of check_path_freq_vs_pi_old

int
path_array_insert(Path* path_array[], int count_array[], Path* path) // like path tree but do it dumb way with array
{
    int i=1;
    int paths_equiv = UNKNOWN;
    int count;
    static int long_count = 0;

    if(path->length > MAXSHORTPATHLENGTH){
        long_count++;
        count_array[0] = long_count;
        return 0;
    }
    do{
        if(path_array[i] == NULL ){
            path_array[i] = copy_path(path);
            count_array[i] = 1;
            printf("adding new short path to array. i = %i %p %i\n", i, path_array[i], path_array[i]->length);
                //  print_path(stdout, path); printf("\n");
            return i;
        }
        else{
            paths_equiv = (compare_paths(path_array[i], path, &count) == 0);
            if(paths_equiv){
                count_array[i] += 1; // increment count
                return i;
            }
            else{
                i++;
            }
        }
    }while(TRUE);
} // end of path_array_insert


inline int
update_lambdaI(const Path* path, Lambdas* L, double width, double lambda_max, double exponent)
{
        // if ITselect is 0 updates lambdaI, else updates lambdaT
        // ITavg should be Iavg if updating lambdaI, else Tavg.
        // returns 1 if update accepted, 0 otherwise
    Lambdas L_prop = *L;
    double log_pi_ratio, log_p_accept;
    int NpM = path->start->perm->n_edge;
    double L_g = path->start->perm->L_g;
    int result; 

        // do a MH update of lambdaI
           
    if(width >= 1.99*lambda_max) width = 1.99*lambda_max; 
    L_prop.lambdaI = fabs(L->lambdaI + width*(drand(&rng_long) - 0.5)); // get prop by adding random bit and reflecting from 0 and lambda_max
    if(L_prop.lambdaI > lambda_max) L_prop.lambdaI = 2.0*lambda_max - L_prop.lambdaI; // reflect from lambda_max
        
    L_prop.Lambda = get_Lambda(L_g, L_prop.lambdaI, L->lambdaT);

  
    log_pi_ratio = (L->Lambda - L_prop.Lambda)
        + log(L_prop.Lambda/L->Lambda)*(double)(path->length_a - path->length)
        + log(L_prop.lambdaI/L->lambdaI)*(double)path->length_i
        + log_prod_a_to_n(path, &L_prop) - log_prod_a_to_n(path, L); // pi_prop/pi_old

    log_p_accept = log_pi_ratio*exponent;
    if(!is_normal_or_zero(log_p_accept)) n_lambdaIT_p_a_nnoz++;
    if((log_p_accept > 0.0) || log(drand(&rng_long)) < log_p_accept){ // accept proposal
        L->lambdaI = L_prop.lambdaI;
        L->Lambda = L_prop.Lambda;
            // L->lambdaT is unchanged
        L->xi = get_xi(L->lambdaI, L->lambdaT);
        L->r = (L->lambdaI/L->lambdaT);
        result = 1;
    }
    else result = 0; // proposal rejected
       
    return result;
}

inline int
update_lambdaT(const Path* path, Lambdas* L, double width, double lambda_max, double exponent)
{
        // if ITselect is 0 updates lambdaI, else updates lambdaT
        // ITavg should be Iavg if updating lambdaI, else Tavg.
        // returns 1 if update accepted, 0 otherwise
    Lambdas L_prop = *L;
    double log_pi_ratio, log_p_accept;
    int NpM = path->start->perm->n_mark + path->start->perm->n_chrom;
    double L_g = path->start->perm->L_g;
    int result;

// do a MH update of lambdaT
              
    if(width >= 1.99*lambda_max) width = 1.99*lambda_max; 
    L_prop.lambdaT = fabs(L->lambdaT + width*(drand(&rng_long) - 0.5)); // get prop by adding random bit and reflecting from 0 and lambda_max
    if(L_prop.lambdaT > lambda_max) L_prop.lambdaT = 2.0*lambda_max - L_prop.lambdaT; // reflect from lambda_max
        
    L_prop.Lambda = get_Lambda(L_g,  L->lambdaI, L_prop.lambdaT);
   
    log_pi_ratio = (L->Lambda - L_prop.Lambda)
        + log(L_prop.Lambda/L->Lambda)*(double)(path->length_a - path->length)
        + log(L_prop.lambdaT/L->lambdaT)*(double)path->length_t
            //  + log(prod_a_to_n(path, &L_prop)) - log(prod_a_to_n(path, L)); // pi_prop/pi_old
      + log_prod_a_to_n(path, &L_prop) - log_prod_a_to_n(path, L); // pi_prop/pi_old   
   
    log_p_accept = log_pi_ratio*exponent; //
    if(!is_normal_or_zero(log_p_accept)) n_lambdaIT_p_a_nnoz++;
    if((log_p_accept > 0.0) || log(drand(&rng_long)) < log_p_accept){ // accept proposal
            // L->lambdaI unchanged
        L->Lambda = L_prop.Lambda;
        L->lambdaT = L_prop.lambdaT;
        L->xi = get_xi(L->lambdaI, L->lambdaT);
        L->r = (L->lambdaI/L->lambdaT);
        result = 1;
    }
    else result = 0; // proposal rejected
    
    return result;
}


inline int
update_r(Path* path, Lambdas* L)
{
        //  double pi_old, pi_prop, p_accept;
    double log_pi_old, log_pi_prop, log_p_accept;
    double r = L->r;
    double log_prior_factor, log_prior_factor_p;
    Lambdas L_prop = *L;
    int result;
    double L_g = path->start->perm->L_g;

    // fprintf(stderr, "top of update_r. r: %g\n", r);
        // MH step for r
                   
    L_prop.r = fabs(r + r_STEP_WIDTH*(drand(&rng_long) - 0.5));  // reflect about r = 0. For this proposal dist, q_ratio is 1.0
    L_prop.xi = (L_prop.r <= 1.0)? L_prop.r: 2.0 - 1.0/L_prop.r;
    log_prior_factor = (r < 1.0)? 0.0: -2.0*log(r); // from prior for r
    log_prior_factor_p = (L_prop.r < 1.0)? 0.0: -2.0*log(L_prop.r); // from prior for r

    //fprintf(stderr, "QQQQQ: %g  %g  %g  %g \n", L_prop.r, L_prop.xi, log_prior_factor, log_prior_factor_p);

    get_lambdaIT_from_Lambda_r(L_g, L->Lambda, L_prop.r, &(L_prop.lambdaI), &(L_prop.lambdaT));
    
    log_p_accept = log(L_prop.lambdaI/L->lambdaI)*(double)path->length_i + log(L_prop.lambdaT/L->lambdaT)*(double)path->length_t;
    //printf("A %7.5f\n", log_p_accept);
    log_p_accept += log_prod_a_to_n(path, &L_prop);
    //printf("B %7.5f\n", log_p_accept);

    log_p_accept += log_prior_factor_p;
    //printf("C %7.5f\n", log_p_accept);

    log_p_accept -= log_prod_a_to_n(path, L);
    //printf("D %7.5f\n", log_p_accept);

    log_p_accept -= log_prior_factor;
    //printf("E %7.5f\n", log_p_accept);


    if(!is_normal_or_zero(log_p_accept)){
   /*      printf("L, L_prop: %g %g %g %g %g  %g %g %g %g %g\n", L->Lambda, L->lambdaI, L->lambdaT, L->r, L->xi, */
/*                L_prop.Lambda, L_prop.lambdaI, L_prop.lambdaT, L_prop.r, L_prop.xi);    */
        printf("in update_r. log_p_accept is not normal. log_p_accept= %g \n", log_p_accept); n_rxi_p_a_nnoz++;}
    
    if(log_p_accept > 0.0 || log(drand(&rng_long)) < log_p_accept){ // accept proposal     
        *L = L_prop;
        result = 1;
    }
    else result = 0;
    return result;
}


inline int
update_xi(Path* path, Lambdas* L)
{
    double log_pi_old, log_pi_prop, log_p_accept;
    double r = L->r, xi = L->xi;
    Lambdas L_prop = *L;
    int result;
    double L_g = path->start->perm->L_g;
       
        // use xi
    assert(XI_WIDTH < 4.0);
    L_prop.xi = fabs(xi + XI_WIDTH*(drand(&rng_long) - 0.5));
    L_prop.xi = ( L_prop.xi >= 2.0)? 4.0 -  L_prop.xi: L_prop.xi; // now L_prop.xi should be between 0 and 2
    L_prop.r = ( L_prop.xi <= 1.0)?  L_prop.xi: 1.0/(2.0 -  L_prop.xi);

    get_lambdaIT_from_Lambda_r(L_g, L->Lambda, L_prop.r, &(L_prop.lambdaI), &(L_prop.lambdaT));
    
    log_p_accept = log(L_prop.lambdaI/L->lambdaI)*(double)path->length_i + log(L_prop.lambdaT/L->lambdaT)*(double)path->length_t;
    log_p_accept += log_prod_a_to_n(path, &L_prop);
    log_p_accept -= log_prod_a_to_n(path, L);

    if(!is_normal_or_zero(log_p_accept)){ printf("in update_xi. log_p_accept is not normal. p_accept= %g \n", log_p_accept);  n_rxi_p_a_nnoz++;}   
      
    if(log_p_accept > 0.0 || log(drand(&rng_long)) < log_p_accept){ // accept proposal    
        *L = L_prop;
        result = 1;
    }
    else result = 0; // reject proposal
    return result;
}


// end of chain.c
