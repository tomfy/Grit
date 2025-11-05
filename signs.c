// signs.c functions for dealing with signs of markers in unsigned case
#include "grit.h"
#include "structs.h"
#include "perm.h"
#include "exts.h"
#include "signs.h"
#include "rngs.h"
#include "gr_misc.h"
#include "params.h"


double
flip_signs(const Run_info_in* r_in, Permutation* perm1old, Permutation* perm1new, Permutation* perm2, int* flips)
{
        // flips signs in perm1new and returns probability ratio for flips q_signs_backward/q_signs_forward
    int n;
    double flip_probs[NMARKERMAX+1]; 
    double q_signs_forward, q_signs_backward, q_signs_ratio = -1;
    int choose_signs = r_in->Choose_signs;
        
    if(choose_signs == FLIPUSECDMIX){
        if(drand(&rng_long) < 0.5) { 
            choose_signs = FLIPONEUSECD;
        }
        else {
            choose_signs = FLIPSOMEUSECD;
        }
    }

    if(choose_signs == FLIPSOMEUSECD){
        printf("in flip_signs. option choose_signs = %i not implemented. Setting choose_signs to %i . \n", FLIPSOMEUSECD, FLIPONEUSECD);
        choose_signs = FLIPONEUSECD; getchar();
    }
    if(choose_signs == FLIPONEUSECD){
        q_signs_forward = choose_one_to_flip(perm1old, perm2, r_in->Flip_one_sign_epsilon, flips); // choose 1 flip, return probability. 
        assert(real_close(q_signs_forward, prob_one_to_flip(perm1old, perm2, r_in->Flip_one_sign_epsilon, flips), DBL_EPSILON*RCFACTOR_1));  
        flip_signs_by_position(perm1new, flips); // do length 1 inversions to get perm1new
            // printf("in flip_signs, after flip_signs_by_position\n");

        q_signs_backward = prob_one_to_flip(perm1new, perm2, r_in->Flip_one_sign_epsilon, flips); // prob of flip, backwards
            //   printf("q_signs, back, forw: %g %g \n", q_signs_backward,q_signs_forward);
        q_signs_ratio = q_signs_backward/q_signs_forward; 
    }
    else if(choose_signs == FLIPSOMEUSECD){
        printf("in flip_signs. option choose_signs = %i not implemented.  (try choose_signs = %i ). \n", FLIPSOMEUSECD, FLIPONEUSECD);
      /*   n = perm1new->n_mark + perm1new->n_chrom - 1;    */
/*         get_flip_probs(perm1old, perm2, r_in->Flip_some_sign_epsilon, DELTADM1_FLIP_PROB,  flip_probs);            */
/*         q_signs_forward = get_flips(n, flip_probs, flips); // throw the dice, decide which to flip */
/*         flip_signs_by_position(perm1new, flips); // actually do length 1 inversions to get new signs for perm1new_ */
/*         assert(real_close(q_signs_forward, get_flips_prob(n, flip_probs, flips))); */
/*         get_flip_probs(perm1new, perm2, r_in->Flip_some_sign_epsilon, DELTADM1_FLIP_PROB, flip_probs); */
/*         q_signs_backward = get_flips_prob(n, flip_probs, flips); */
/*         q_signs_ratio = q_signs_backward/q_signs_forward;  */
    }
    else if(choose_signs == FLIPONERANDOM){ 
        flip_one(perm1new, flips); // just flip one of the markers
	   printf("in flip_signs. option choose_signs = %i not implemented.  (try choose_signs = %i ). \n", FLIPONERANDOM, FLIPONEUSECD);
    }
    else{
        printf("unknown Choose_signs option: %i \n", choose_signs);
        process_error();
    }
    if(q_signs_ratio < 0) exit(1);
    return q_signs_ratio;
} // end of flip_signs


void
get_flip_deltads(Permutation* perm1, Permutation* perm2, int* deltads, int* seds)
{ // get cycle decomposition, set deltad[i] = the change in d (delta_d = delta_c_g - delta_c_i) that would result
        // from flipping the marker in position i+1
    int i;
    Cycle_decomposition the_cd;
    Cycle_element* edge1, * edge2;
    Cycle* eds1, *eds2;
   
    get_CD(perm1, perm2, &the_cd);
        //   print_edss(&the_cd);
    for(i=0; i<perm1->n_mark + perm1->n_chrom - 1; i++){    //
        edge1 = the_cd.elements + i;
        edge2 = edge1 + 1;
        assert(i == edge1->position);
        assert(edge2->position == edge1->position+1);
        if(edge1->chromosome == edge2->chromosome){ // this pair of edges are the neighbors of a real (internal) marker
            if(edge1->eds == edge2->eds){ // same eds, d_c_g = 0, d_c_i = 0 or 1
                seds[i] = TRUE;
                if(edge1->direction == edge2->direction){ deltads[i] = 0;}
                else{ deltads[i] = -1;}
            }
            else{  // distinct edss
                seds[i] = FALSE;
                eds1 = the_cd.edss + edge1->eds; eds2 = the_cd.edss + edge2->eds;
                if(eds1->is_end_cycle && eds2->is_end_cycle){ // both edss have ends. delta c_i = 0
                    if((eds1->first_black + eds1->last_black + eds2->first_black + eds2->last_black) != 2){ // bg+bb, bb+bb, etc. delta c_g = 0
                        deltads[i] = 0;  //0
                    }
                    else{ // delta c_g is -1, 0 or 1
                        if(eds1->first_black == eds1->last_black){ // bb + gg, or gg + bb, delta c_g = -1
                            deltads[i] = -1;
                        }
                        else{ // 
                            if((eds1->first_black == eds2->first_black) == (edge1->direction == edge2->direction)){
                                    // bg+bg, same dir, bg+gb opp. dir., etc. delta c_g = 1
                                deltads[i] = 1;
                            }
                            else{  // bg+bg opp. dir, bg+gb same dir, etc. delta c_g = 0
                                deltads[i] = 0;
                            }
                        }
                    }  
                }
                else{ // at least one edge is on int cycle -> d_c_i = -1, d_c_g = 0;
                    deltads[i] = 1; // step away
                }
            }
        }
        else{
            seds[i] = FALSE; // just a swap of chromosome end markers - doesn't help with hurdles
            assert(edge1->at_right_end && edge2->at_left_end);
            deltads[i] = NOTFLIP;
        }
            //  printf("i, deltads[i]: %i %i \n", i, deltads[i]);
    }
} // end of get_flip_deltads


void
get_flip_probs(Permutation* p1, Permutation* p2, double sign_epsilon, double deltadm1_flip_prob, double* flip_prob)
{
        // for FLIPSOMEUSECD
        // deltads[i] is delta d that would result from flipping
        // marker in position i+1. n is the number of markers
    int i, the_delta_d;
        //double sign_epsilon = r_in->Sign_epsilon;
    int deltads[NMARKERMAX+1];
    int seds[NMARKERMAX+1];

    get_flip_deltads(p1, p2, deltads, seds);
    
    for(i=0; i<p1->n_mark + p1->n_chrom - 1; i++){
        the_delta_d = deltads[i];
        if(the_delta_d == 1) { // away
                //flip_prob[i] = sign_epsilon/(1.0 + sign_epsilon);
                flip_prob[i] = sign_epsilon;
        } 
        else if(the_delta_d == 0) { // delta d = 0
            flip_prob[i] = DELTAD0_FLIP_PROB;
        }
        else if(the_delta_d == -1) { // toward target
                //flip_prob[i] = 1.0/(1.0 + sign_epsilon);
                 flip_prob[i] = deltadm1_flip_prob;
        }
        else {
            printf("in get_flip_probs. delta_d has bad value: %i \n", the_delta_d);
            process_error();
        } 
    } 
} // end of get_flip_probs

 
double 
get_flips(int n, double* flip_probs, int* flips) 
{ 
        // throw the dice, set flips array according to probs in flip_probs 
        // returns the probability of this set of flips 
        // flips[i] = 1  ->  marker at position i+1 is to be flipped 
     
    int i; 
    double the_flip_prob, prob = 1.0; 
     
    for(i = 0; i < n; i++){ 
        the_flip_prob = flip_probs[i]; 
        if(drand(&rng_long) < the_flip_prob) { 
            flips[i] = 1; 
            prob *= the_flip_prob; 
        } 
        else { 
            flips[i] = 0; 
            prob *= 1.0 - the_flip_prob; 
        }
    } 
    return prob; 
} // end of get_flips
 
void 
fix_signs(Path* path, Step* last, int* flips_by_number)  
{ 
        // run along path up to but not including last, flipping signs 
        //  
    Step* the_step = path->start; 
     
    do{ 
        if(the_step == last) { printf("In fix_signs; the_step == last\n"); process_error(); }
        assert(the_step != last); 
        flip_signs_by_number(the_step->perm, flips_by_number);
            // also need to fix the rev's of the steps
        fix_rev(the_step->perm, &(the_step->rev), flips_by_number);  
        the_step = the_step->next; 
    }while(the_step != last); 
}         

void
fix_rev(Permutation* perm, Reversal* rev, int* flips_by_number)
{
    int i, mnumb;
    int* numbs;
    int* marker_numbers;
    int checks;

 /*    printf("top of fix rev \n"); */
/*     printf("rev->numbs: %i %i %i %i \n", */
/*            rev->numbs[0], rev->numbs[1], */
/*            rev->numbs[2], rev->numbs[3]); */
/*     printf("rev->marker_numbers: %i %i %i %i \n", */
/*            rev->marker_numbers[0], rev->marker_numbers[1], */
/*            rev->marker_numbers[2], rev->marker_numbers[3]); */
        //  check_rev(16, rev);
    for(i=0, numbs = rev->numbs; i<4; i++, numbs++){
        if(*numbs != UNKNOWN){ // skip end markers (numbs[i] UNKNOWN)
            mnumb = (*numbs +1)/2; // get marker number corresponding to marker end number *numb
            if(flips_by_number[mnumb-1]) {
                *numbs = 4*mnumb - 1 - *numbs; // set *numbs to number of other marker end of same marker
            }
        }
    }

    for(i=0, numbs = rev->marker_numbers; i<4; i++, numbs++){
        if(!perm->p[*numbs].is_end_marker){ // skip end markers 
            mnumb = (*numbs +1)/2; // get marker number corresponding to marker end number *numb
            if(flips_by_number[mnumb-1]) {
                *numbs = 4*mnumb - 1 - *numbs; // set *numbs to number of other marker end of same marker
            }
        }
    }
    
  /*   checks = check_rev(perm->n_mark, rev); */
/*         //  printf("bottom of fix rev \n"); */
/*     if(!checks) getchar(); */
    assert( check_rev(perm->n_mark, rev));
}

double 
get_flips_prob(int n, double* flip_probs, int* flips) 
{ 
        // get prob of flips given flip_probs 
        // flips is an array of 1 (-> flip) and 0 (-> don't flip) elements 
        // flip_probs an array of the probs of flipping each marker 
    int i; 
    double prob = 1.0;
     
    for(i=0; i<n; i++){ 
        if (flips[i] == TRUE) prob *= flip_probs[i]; 
        else prob *= 1.0 - flip_probs[i]; 
    } 
    return prob; 
} // end of get_flips_prob


double
choose_one_to_flip(Permutation* perm1, Permutation* perm2, double sign_epsilon, int* flips)
{
        // deltad[i] is the delta d (-1,0,+1) that would result from flip of marker in position i+1
        // or if deltad[i] = NOTFLIP then position i is not a marker (it is a pair of marker ends on neighboring chromosomes)
        // choose one of these 3 possibilities (1,0,-1), with rel probs eps4^2, eps4, 1, respectively
        // choose one of the flips of this kind uniformly at random, store in flips array
        // return probability of this choice of marker to flip

    int i, nplus=0, nzero=0, nzero_seds=0, nzero_deds=0, nminus=0, nnotaflip=0, the_delta_d, ipick;
    int plusoneflips[NMARKERMAX+MAX_N_CHROMOSOMES], zeroflips[NMARKERMAX+MAX_N_CHROMOSOMES], minusoneflips[NMARKERMAX+MAX_N_CHROMOSOMES];
    int zeroflips_seds[NMARKERMAX+MAX_N_CHROMOSOMES], zeroflips_deds[NMARKERMAX+MAX_N_CHROMOSOMES]; 
    double p1, p0, pm1, sump;
        //double sign_epsilon = r_in->Sign_epsilon;
    double this_flip_prob, the_rand_number;
    int deltads[NMARKERMAX+MAX_N_CHROMOSOMES];
    int seds[NMARKERMAX+MAX_N_CHROMOSOMES];
    double pdd0_seds;

    get_flip_deltads(perm1, perm2, deltads, seds);
    
    for (i=0; i<perm1->n_mark + perm1->n_chrom - 1; i++){
        flips[i] = FALSE; // initialize flips array
        the_delta_d = deltads[i];
        if (the_delta_d == -1){
            minusoneflips[nminus] = i;
            nminus++;
        }
        else if (the_delta_d == 0){
            zeroflips[nzero] = i;
            nzero++;
            if(seds[i] == TRUE){
                zeroflips_seds[nzero_seds] = i;
                nzero_seds++;
            }
            else{
                zeroflips_deds[nzero_deds] = i;
                nzero_deds++;
            }
        }
        else if (the_delta_d == 1){
            plusoneflips[nplus] = i;
            nplus++;
        }
         else if (the_delta_d == NOTFLIP){
            nnotaflip++;
        }
        else {
            printf("in choose_one_to_flip, the_delta_d has bad value: %i \n", the_delta_d);
            process_error();
        }
    }

    assert(nnotaflip == perm1->n_chrom - 1); 
    
        //nzerop = nzero + ADD_TO_NZERO; // add 1 or 0
// set probs to zero if no flips give that deltad
    p1 = GAMMA*sign_epsilon*sign_epsilon*step_fcn(nplus);
        //  printf("in choose_one_to_flip. GAMMA: %g \n", GAMMA);
    p0 = sign_epsilon*(step_fcn(nzero));
    pm1 = step_fcn(nminus); 
    sump = p1 + p0 + pm1; // normalize
    p1 /= sump; p0 /= sump; pm1 /= sump;

        // throw the dice
    the_rand_number = drand(&rng_long); // just get rand number once; don't call drand in each branch!!!!!
    if(the_rand_number < pm1){
            // choose a deltad=-1 flip
        ipick = (int)(drand(&rng_long)*nminus); // in range 0 to nminus-1
        flips[minusoneflips[ipick]] = TRUE; 
        this_flip_prob = pm1/(double)nminus;
    }
    else if (the_rand_number < pm1 + p0){
            // choose a deltad=0 flip, or no flip (also deltad=0)
        if((nminus > 0 ) || !PREFERDD0SEDS_SIGNFLIP){
            ipick = (int)(drand(&rng_long)*(nzero)); // in range 0 to nzero - 1 -> flip; nzero -> no flip
            if (ipick < nzero) { flips[zeroflips[ipick]] = TRUE;}
                // else ipick == nzero don't set an elem of flips to TRUE; -> no flip
            this_flip_prob = p0/(double)(nzero);
        }
        else{ // prefer same eds delta_d=0 flips
                //  printf("nzeros: %i %i %i \n", nzero, nzero_seds, nzero_deds);
           /*  if(nzero_seds >= nzero_deds) */
/*                 printf("nzero, nzero_seds, nzero_deds : %i %i %i **** \n", nzero, nzero_seds, nzero_deds); */
/*             else  printf("nzero, nzero_seds, nzero_deds : %i %i %i \n", nzero, nzero_seds, nzero_deds); */
            if(nzero_deds == 0){
                assert(nzero_seds > 0);
                pdd0_seds = 1.0;
            }
            else if(nzero_seds == 0){
                assert(nzero_deds > 0);
                pdd0_seds = 0.0;
            }
            else{
                assert(nzero_seds > 0 && nzero_deds > 0);
                pdd0_seds = PDD0SEDS_SIGNFLIP;
            }
                
            if(drand(&rng_long) < pdd0_seds){ // do a seds delta_d=0 flip
                    //  printf("In choose_one_to_flip. seds branch\n");
                ipick = (int)(drand(&rng_long)*(nzero_seds)); // 
                if (ipick < nzero_seds) { flips[zeroflips_seds[ipick]] = TRUE;}
                    // else ipick == nzero don't set an elem of flips to TRUE; -> no flip
                this_flip_prob = pdd0_seds*p0/(double)(nzero_seds);
            }
            else{ // do a distinct eds delta_d=0 flip
                    //   printf("In choose_one_to_flip. deds branch\n");
                ipick = (int)(drand(&rng_long)*(nzero_deds)); // 
                if (ipick < nzero_deds) { flips[zeroflips_deds[ipick]] = TRUE;}
                    // else ipick == nzero don't set an elem of flips to TRUE; -> no flip
                this_flip_prob = (1.0-pdd0_seds)*p0/(double)(nzero_deds);
            }
        }
    }
    else{
            // choose a deltad=1 flip
        ipick = (int)(drand(&rng_long)*nplus); // in range 0 to nplus-1
        flips[plusoneflips[ipick]] = TRUE;
        this_flip_prob = p1/(double)nplus;
    }
   
    return this_flip_prob;
} // end of choose_one_to_flip


double
prob_one_to_flip(Permutation* perm1, Permutation* perm2, double sign_epsilon, int* flips)
{
        // given flips (flip of one marker), get deltads and get probability.
        // deltads[i] is delta d due to flip of marker in position i+1
        // flips[i] == TRUE -> flip marker at position i+1
        // return probability of this choice of marker to flip

    int i, nplus=0, nzero=0, nzero_seds=0, nzero_deds=0, nminus=0, nnotaflip = 0, the_delta_d, nflips=0, the_one_to_flip=-1;
    double p1, p0, pm1, sump;
    double this_flip_prob;
    int deltads[NMARKERMAX+MAX_N_CHROMOSOMES];
    int seds[NMARKERMAX+MAX_N_CHROMOSOMES];

    get_flip_deltads(perm1, perm2, deltads, seds);

    for (i=0; i<perm1->n_mark + perm1->n_chrom - 1; i++){
        if(flips[i]){
            nflips++;
            the_one_to_flip = i; // 0 to N-1
        }
        the_delta_d = deltads[i];
        if (the_delta_d == -1){ nminus++; }
        else if (the_delta_d == 0){
            nzero++;
            if(seds[i] == TRUE) nzero_seds++;
            else nzero_deds++;
        }
        else if (the_delta_d == 1){ nplus++; }
        else if (the_delta_d == NOTFLIP){ nnotaflip++; }
        else { printf("in prob_one_to_flip, the_delta_d has bad value: %i \n", the_delta_d); process_error();}
    }
    if(nflips == 1){   
        the_delta_d = deltads[the_one_to_flip];
        assert(the_delta_d != NOTFLIP);
    }
    else if(nflips == 0){
        printf("in prob_one_to_flip. no flip done. \n");
        getchar();
        the_delta_d = 0;
    }
    else {
        printf("in prob_one_to_flip. nflips has bad value: %i \n", nflips);
        process_error();
    }
        // set probs to zero if no flips give that deltac
        // except "no flip" is an option (delta c = 0) so p0 never zero
        // nzerop = nzero + ADD_TO_NZERO;
    p1 = GAMMA*sign_epsilon*sign_epsilon*step_fcn(nplus); // delta d = 1 is suppressed by ~ sign_epsilon^2
    p0 = sign_epsilon*(step_fcn(nzero)); // delta d = 0 is suppressed by ~ sign_epsilon
    pm1 = step_fcn(nminus); // delta d = -1 is preferred rel to delta d of 1 or 0
    sump = p1 + p0 + pm1; // normalize
    p1 /= sump; p0 /= sump; pm1 /= sump;
       
    if(the_delta_d == 1){
        this_flip_prob = p1/(double)nplus;
    }
    else if (the_delta_d == 0){
        if((nminus > 0) || !PREFERDD0SEDS_SIGNFLIP){
        this_flip_prob = p0/(double)(nzero);
        }
        else{ // prefer seds
            if(nzero_seds == 0){
                assert(nzero_deds > 0);
                assert(seds[the_one_to_flip] == FALSE);
                this_flip_prob = p0/(double)nzero_deds;
            }
            else if(nzero_deds == 0){
             assert(nzero_seds > 0);
                assert(seds[the_one_to_flip] == TRUE);
                this_flip_prob = p0/(double)nzero_seds;
            }
            else{
                if(seds[the_one_to_flip] == TRUE){
                    this_flip_prob = PDD0SEDS_SIGNFLIP*p0/(double)nzero_seds;
                }
                else{
                    this_flip_prob = (1.0 - PDD0SEDS_SIGNFLIP)*p0/(double)nzero_deds;
                }
            }
        }           
    }
    else {
        this_flip_prob = pm1/(double)nminus;
    }
    return this_flip_prob;
} // end of prob_one_to_flip


double 
randomize_signs(int n, double pflip, int* flipped) 
{ 
        // set each elements of flipped to TRUE with probability  pflip 
        // returns probability of signs chosen 
    int i; 
    double prob = 1.0; 
     
    for(i = 0; i < n; i++){ 
        if(drand(&rng_long) < pflip) {
            flipped[i] = TRUE; prob *= pflip;
        }
        else {
            flipped[i] = FALSE; prob *= (1.0 - pflip);
        }
    } 
    return prob; 
} 
 
void 
flip_one(Permutation* perm, int* flips) 
{ 
        // just flip one marker at random 
    int i, n = perm->n_mark + perm->n_chrom, l_cut;
    Reversal rev;
    
    l_cut = (int)(drand(&rng_long)*(n-1));
    rev = make_signflip_rev(perm, l_cut, l_cut+1);
    reverse(perm, &rev);
    printf("In flip_one. after reverse\n");
    for(i=0; i<n; i++){ 
        flips[i] = FALSE; 
    } 
    flips[l_cut] = TRUE;
} 
 
void 
flip_signs_by_position(Permutation* perm, int* flips) 
{ 
        // flip signs of markers in perm as specified by flips array 
        // which indicates POSITIONS of markers to be flipped 
        // i.e. flips[3] = 1 -> flip marker 4 (cutting black edges 3 and 4) 
         
     
    int n = perm->n_mark + perm->n_chrom, l_cut;  // this is the number of the left black edge to be cut to flip the marker 
    Reversal rev;
    
    for(l_cut = 0; l_cut < n-1; l_cut++){ 
        if(flips[l_cut]) {
            rev = make_signflip_rev(perm, l_cut, l_cut+1);
            reverse(perm,  &rev);
            // printf("In flip_signs by position. after reverse\n");
        }
    }
            //  printf("before returning from flip signs by poistion\n");
} 


void 
flip_signs_by_number(Permutation* perm, int* flips) 
{ 
        // flips is here an array specifying whether the marker with a certain number gets flipped 
        // specifically flips[i] = 1 -> flip marker number i+1 
    int i, mpos, n = perm->n_mark + perm->n_chrom; // 
    Marker_end* m;
    Reversal rev;

    n = perm->n_mark; // just flip real markers 
 
    for(i=0; i<n; i++){
       m = perm->p + (2*i+1); // do marker_ends 1, 3, ... 
       mpos = (m->position + 1)/2; 
       if (flips[i]) {
           rev = make_signflip_rev(perm, mpos-1, mpos);
           reverse(perm, &rev);
               //  printf("In flip_signs by number. after reverse\n");
       }
    }
}


void 
get_flips_by_number(Permutation* perm, int* flips_by_pos, int* flips_by_number) 
{ 
        // flips_by_pos[i] = 1 means marker at position i+1 is to be flipped 
        // flips_by_number[i] = 1 means marker with number i+1 is to be flipped 
    int n = perm->n_mark + perm->n_chrom - 1; // this is the number of edges 
    int i, mpos, flip, mnumber; 
    Marker_end* m = perm->p; 
 
    for (i=0, m = m->black; i<n; i++, m = m->other->black){ 
        mpos = (m->position + 1)/2; // marker position (1 through N) 
        if(i != (mpos - 1)) { printf("prob in get_flips_by_number. i, mpos-1 (should be equal): %i %i \n", i, mpos-1); process_error();}
        flip = flips_by_pos[mpos - 1]; // whether it should be flipped or not 
        mnumber = (m->markend_numb + 1)/2; // marker number (1 through N) 
        flips_by_number[mnumber - 1] = flip; // whether it should be flipped or not 
    } 
}


int get_signs_short(Permutation* p1, Permutation* p2, const Run_info_in* r_in)
{
        // flips signs of markers in p1
        // so as to shorten distance from p1 to p2
    int iii, iiii, flips[NMARKERMAX];
    double factor;
    int n_T = 5; // number of times through iiii loop
    Cycle_decomposition the_cd;
        //   double flip_probs[NMARKERMAX];
   
 printf("in get_signs_short. top. \n");
    (void)randomize_signs(p1->n_mark, 0.5, flips);
        //  flip_signs_by_position(p1, flips);
    flip_signs_by_number(p1, flips);
printf("in get_signs_short. after first flip_signs_by_number call. \n");
        // fix signs of conserved segments; singletons signs left unchanged (i.e. random)
           get_conserved_flips(p1, p2, flips);
        printf("in get_signs_short. after get_conserved_flips. \n");
           flip_signs_by_number(p1, flips);
                
 /*    get_flip_probs(p1, p2, r_in->Flip_some_sign_epsilon, 1.0 - r_in->Flip_some_sign_epsilon, flip_probs);
// this to get good flips all at once (i.e. FLIPSOMEUSECD style)            */
           printf("in get_signs_short. after flip .. by number. \n");                
/*     (void)get_flips(p1->n_mark, flip_probs, flips); // throw the dice, decide which to flip */
/*     flip_signs_by_position(p1, flips); */

    for(iiii=0; iiii<n_T; iiii++){ // something like sim annealing; turn down the T each time through loop
        factor = (n_T - iiii)/(double)n_T;
        for(iii=0; iii<2*p1->n_mark; iii++) // this loop to get good flips one at a time (i.e. FLIPONEUSECD style)
        {
            (void)choose_one_to_flip(p1, p2, factor*r_in->Flip_one_sign_epsilon, flips); // choose 1 flip, return probability.
            flip_signs_by_position(p1, flips);
                //  printf("iiii, iii, in get_signs_short, after flip_signs_by_position. %i %i \n", iiii, iii);
        }
    }
    get_CD(p1, p2, &the_cd);
    return the_cd.n_mark - the_cd.n_chrom - the_cd.n_int_cycles + the_cd.c_g; // this is the "distance" from p1 to p2 
} // end of get_signs_short

int
get_conserved_flips(Permutation* perm1, Permutation* perm2, int* flips)
{
        // look at perm1, perm2, store in flips the marker numbers of flips necessary
        // to produce most likely relative signs for conserved segments
        //

    Marker_end* end1, * end2;
    int i, m1, m2;
    int pos1, pos2;
    int n_flips = 0;

    for(i = 0; i<perm1->n_mark; i++){
        flips[i] = FALSE;
    }
    for(i=0, end1 = perm1->p; end1 != NULL; i++, end1 = end2->other){
        end2 = end1->black;
        assert(end2->position == end1->position+1);
        assert(end1->position % 2 == 0);
        m1 = end1->markend_numb; m2 = end2->markend_numb;
        pos1 = perm2->p[m1].position; pos2 = perm2->p[m2].position;
        if(!(end1->is_end_marker || end2->is_end_marker)){
            if((pos2+1)/2 == (pos1+1)/2+1){ // pair of markers is in same order in both genomes, signs should be same too
                if(pos1 % 2 != 0) { flips[(m1-1)/2] = TRUE; n_flips++;}
                if(pos2 % 2 == 0) { flips[(m2-1)/2] = TRUE; n_flips++;}
            }
            else if((pos1+1)/2 == (pos2+1)/2+1){
                if(pos1 % 2 == 0) { flips[(m1-1)/2] = TRUE; n_flips++;}
                if(pos2 % 2 != 0) { flips[(m2-1)/2] = TRUE; n_flips++;}
            }
        }
            // else markers are not adjacent in perm2 - do nothing
    }
    return n_flips;
} // end of get_conserved_flips
            
            

// end of signs.c
