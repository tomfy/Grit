// path_etc.c  functions for allocation, checking, printing, etc.
// related to paths

#include "grit.h"
#include "params.h"
#include "structs.h" 
#include "perm.h"
#include "perm_etc.h"
#include "exts.h"
#include "path.h"
#include "path_etc.h"
#include "gr_misc.h"
#include "chain.h"

//#include "path_tree_etc.h"

void freepath(Path* path) 
{ 
   // doesn't free memory assoc with start, end, prob, length. 
  Step* the_step = path->start; 
  Step* the_next_step;
  int count_freed_steps = 0;
 
  while (the_step != path->end) { 
      the_next_step = the_step->next; 
      //assert((the_next_step == NULL) == (the_step == path->end)); 
      permfree(&(the_step->perm)); 
      free(the_step); n_step_alloc--;
            count_freed_steps++;
      the_step = the_next_step; 
  } 
  permfree(&(the_step->perm)); 
  free(the_step); n_step_alloc--; 
  count_freed_steps++;
      //     printf("in freepath. %i steps freed. \n", count_freed_steps); //getchar();
  path->start = NULL; 
  path->end = NULL; 
  path->prob = UNKNOWN; 
  path->length = UNKNOWN;  
} // end of function freepath

Path* path_alloc(Step* step)
{
        // returns a Path* pointing to newly allocated memory
        // which is initialized to be a path of length zero
        // step is a pointer to the first step on path;
        // it doesn't need to be allocated before calling path_alloc
    
    Path* path = (Path*)chmalloc(sizeof(Path)); n_path_alloc++;
    path->accept = UNKNOWN;
    path->start = step; 
    path->end = step;
    path->prob = 1.0; 
    path->length = 0;
    path->length_a = 0;
    path->length_i = 0;
    path->length_t = 0;
    path->length_f = 0;
    path->fission_factor = 1.0;
    path->swap_end_factor = 1.0;
    path->L = 0;
    path->L_i = 0;
    path->L_t = 0;

    return path;
}

Step* step_alloc(Permutation* perm)
{
    Step* step = (Step*)chmalloc(sizeof(Step)); n_step_alloc++; 
    step->perm = perm;
    step->next = NULL;
    step->prev = NULL;
        //   step->rev.left = step->rev.right =
  /*   step->rev.flip_chromosome = step->rev.is_inversion = UNKNOWN;          */
/*     step->rev.numbs[0] = step->rev.numbs[1] = step->rev.numbs[2] = step->rev.numbs[3] = UNKNOWN; */
    step->rev = r_init;
    step->n_dum = 0;

    return step;
}


double
get_fission_factor(const Path* path)
{
        // return fission_factor, 
   
    double fission_factor = 1.0;
    Step* the_step;
    int LLLL = 0;
    
        //    printf("in get_fission_factor\n");
        //    print_path(stdout, path);
    the_step = path->start; 
    while(the_step != path->end){
        if(the_step->rev.is_fission){
            fission_factor *= 2.0*(the_step->perm->n_chrom - the_step->perm->n_chrom_nonempty); // *= 2M_0
                //  printf("fission factor: %g , nchrom, nchrom_ne %i %i \n", fission_factor, the_step->perm->n_chrom, the_step->perm->n_chrom_nonempty);
//  printf("In get_fission_factor. is_fission: %i   ff: %g \n", the_step->rev.is_fission, fission_factor);
        }
            //   printf("LLLL, the_step, the_step->next, path->end: %i %p %p %p \n", LLLL, the_step, the_step->next, path->end);
        the_step = the_step->next; LLLL++;
    }
    return fission_factor; 
} // end of function get_fission_factor

double
get_swap_end_factor(const Path* path)
{
        // return fission_factor, 
   
    double swap_end_factor = 1.0;
    Step* the_step;
    
        //    printf("in get_swap_end_factor\n");
        //    print_path(stdout, path);
    the_step = path->start; 
    while(the_step != path->end){
        if(rev_just_swaps_chromosome_ends(the_step->perm, the_step->rev)){
            swap_end_factor *= pow((double)the_step->perm->n_chrom, 2.0) - (the_step->perm->n_chrom - the_step->perm->n_chrom_nonempty); // *= 2M_0
                //  printf("fission factor: %g , nchrom, nchrom_ne %i %i \n", fission_factor, the_step->perm->n_chrom, the_step->perm->n_chrom_nonempty);
//  printf("In get_fission_factor. is_fission: %i   ff: %g \n", the_step->rev.is_fission, fission_factor);
        }
        the_step = the_step->next;
    }
    return swap_end_factor; 
} // end of function get_swap_end_factor


int
check_path_lengths(const Path* path) 
{ 
        // get path lengths by counting steps until reach path->end
        // check length (length_t, etc) against path->length, etc..
    int length = 0, length_i, length_t = 0, length_f = 0, length_a = 0;
    int result;
    double fission_factor = 1.0;
    Step* the_step;

        //   printf("in check_path_lengths\n");
        //   print_path(stdout, path);
    the_step = path->start; 
    while(the_step != path->end){
        length++;
            //    printf("length, n_dum, length_a: %i %i %i \n", length, the_step->n_dum, length_a);
        length_a += the_step->n_dum; // add # dummy events for this step
        if(!the_step->rev.is_inversion) length_t++;
        if(the_step->rev.is_fission){
            length_f++;
            fission_factor *= 2.0*(the_step->perm->n_chrom - the_step->perm->n_chrom_nonempty); // *= 2M_0
                //  printf("In check_path_lengths: is_fission: %i   ff: %g \n", the_step->rev.is_fission, fission_factor);
        }
        the_step = the_step->next;
    }
    
    length_a += the_step->n_dum; // add in dummies events of last step (after last real inv/trans)
        //  printf("length, n_dum, length_a: %i %i %i \n", length, the_step->n_dum, length_a);
  
    length_a += length;
    length_i = length - length_t;
    if(!(result = ((length == path->length)&&(length_i == path->length_i)
                   &&(length_t == path->length_t)&&(length_f == path->length_f)&&(length_a == path->length_a))))
        printf("in check_path_lengths: \n      length,     length_i, etc: %8i %8i %8i %8i %8i\n  path->length, path->length_i, etc: %8i %8i %8i %8i %8i\n",
               length, length_i, length_t, length_f, length_a, path->length, path->length_i, path->length_t, path->length_f, path->length_a);
    if(!real_close(fission_factor, path->fission_factor, 500.0*DBL_EPSILON)) {
        printf("in check_path_lengths: %g %g \n", fission_factor, path->fission_factor); getchar();
    }
    return result;
}// end of function check_path_lengths


int 
check_path_consistency_multi(const Path* the_path) 
{ 
        // just checks that delta_c_int is -1, 0, or +1 at each step on path
    Step* the_step; 
    Permutation* start_perm, * the_perm, * next_perm; 
    int i, n_inconsistent = 0, delta_c_int;
    Cycle_decomposition the_cd, next_cd;
     
    the_step = the_path->start; 
    start_perm = the_step->perm;

        //   printf("top of check_ath_consistency... \n"); fflush(stdout);
    for(i=0; i<the_path->length; i++){
            //    printf("in check_path_... top of loop. i, %i \n", i); fflush(stdout);
        the_perm = the_step->perm;
            //    printf("the_step->next: %p\n", the_step->next); fflush(stdout);
        next_perm = the_step->next->perm;
        if(!check_rev_perm_consistent(the_perm, &(the_step->rev))) n_inconsistent++;
            //    printf("lskjd \n", i); fflush(stdout);
        get_CD(start_perm, the_perm, &the_cd);
            //    printf("after first get_CD in check_path...\n"); fflush(stdout);
        
        get_CD(start_perm, next_perm, &next_cd);
            //      printf("after second get_CD in check_path...\n"); fflush(stdout);
        delta_c_int = the_cd.n_int_cycles - next_cd.n_int_cycles;
        if(delta_c_int > 1 || delta_c_int < -1){
            print_perm(stdout, start_perm);
            print_perm(stdout, the_perm);
            print_edss(&the_cd);
            print_perm(stdout, next_perm);
            print_edss(&next_cd);
            n_inconsistent++;
        }
            //    printf("in check_ath_consistency...i, the_path->length, the_step, the_step->next %i %i %p %p \n",
            //    i, the_path->length, the_step, the_step->next); fflush(stdout);
        the_step = the_step->next; 
    } 
    if(n_inconsistent > 0) printf("Check_path_consistency_multi found n_inconsistent = %i \n", n_inconsistent); 
    return n_inconsistent; 
} // end of check_path_consistency_multi


void print_dummies(Path* path)
{
        // go along path, print numbers of dummy events for each step
    Step* the_step = path->start, * the_next_step;
    int count = 0, done = FALSE;

    printf("length, length_t, f, a: %i %i %i %i \n", path->length, path->length_t, path->length_f, path->length_a);
 
        // loop over steps on path, including last one
    while(!done) {
        done = (the_step == path->end);
        the_next_step = the_step->next;     assert((the_next_step == NULL) == (the_step == path->end));
        printf("step %i , n_dum: %i \n", count, the_step->n_dum);
        the_step = the_next_step;
        count++;
    }
} // end of function print_dummies

void print_dummies2(Step* step)
{
        // go along path, starting at step, to end printing numbers of dummy events at each step
    Step* the_step = step, * the_next_step;
    int count = 0, done = FALSE;

        //  printf("length, length_t, f, a: %i %i %i %i \n", path->length, path->length_t, path->length_f, path->length_a);
 
        // loop over steps on path, including last one
    while(!done) {
        done = (the_step->next == NULL);
        the_next_step = the_step->next;   //  assert((the_next_step == NULL) == (the_step == path->end));
        printf("step %i , n_dum: %i \n", count, the_step->n_dum);
        the_step = the_next_step;
        count++;
    }
} // end of function print_dummies2

void print_path(FILE* fp, const Path* the_path){
    Step* the_step = the_path->start;
    
    fprintf(fp, "************* printing path: \n");
    while(the_step != the_path->end){
               print_perm(fp,the_step->perm);
        fprintf(fp, "the_step->perm->n_chrom, the_step->perm->n_chrom_nonempty: %i %i \n",
                the_step->perm->n_chrom, the_step->perm->n_chrom_nonempty);
               print_rev(fp, the_step->rev);
        the_step = the_step->next;
    }
           print_perm(fp, the_step->perm);
    fprintf(fp, "the_step->perm->n_chrom, the_step->perm->n_chrom_nonempty: %i %i \n",
                the_step->perm->n_chrom, the_step->perm->n_chrom_nonempty);
    fprintf(fp, "*************** done printing path: \n");
        //  getchar();
}// end of print_path

void print_path_brief(FILE* fp, const Path* the_path)
{
    Step* the_step = the_path->start;
    
    fprintf(fp, "************* printing path: \n");
    while(the_step != the_path->end){
        print_perm_brief(fp,the_step->perm);
        print_rev_brief(fp, the_step->rev); fprintf(fp, "\n");
        the_step = the_step->next;
    }
    print_perm_brief(fp, the_step->perm);

    fprintf(fp, "L, L_i, L_t: %i %i %i   %i %i %i \n", the_path->length, the_path->length_i, the_path->length_t,
            the_path->L, the_path->L_i, the_path->L_t); 
    fprintf(fp, "*************** done printing path: \n");
        //  getchar();
}// end of print_path

void print_dummies1(FILE* fp, const Path* the_path)
{
    Step* the_step = the_path->start;
    
    fprintf(fp, "************* printing dummies on path: \n");
    while(the_step != the_path->end){
        printf("n_dum: %i \n", the_step->n_dum);     
        the_step = the_step->next;
    }
      printf("last step. n_dum: %i \n", the_step->n_dum);         
    fprintf(fp, "*************** done printing dummies on path: \n");
        //  getchar();
}// end of print_dummies


/* void print_path_trans(FILE* fp, Path* the_path) */
/* { // print translocations on path */
/*         // actually just the chromosomes involved in each translocation */
/*     Step* the_step = the_path->start; */
/*     Permutation* the_perm; */
/*     Marker_end* a, * b, * c, * d;  */

/*     fprintf(fp, "%3i %3i  ", the_path->length, the_path->length_t); */
/*     while(the_step != the_path->end){ */
/*         the_perm = the_step->perm; */
/*         get_abcd(the_perm, &(the_step->rev), &a, &b, &c, &d); */
/*         if(a->chromosome != d->chromosome){ */
/*             fprintf(fp, "(%i, %i) ", a->chromosome, d->chromosome); */
/*         }  */
/*         the_step = the_step->next; */
/*     } */
/*     fprintf(fp, "\n"); */
/* }// end of print_path */

void
print_rev(FILE* fp, Reversal rev)
{
    
    fprintf(fp, "inv,toacbd,fiss,seds: %6i %6i %6i %6i  \n dbreakod:%g %g  %g %g   \n marker_numbers:%6i %6i %6i %6i rev.ld,rd:%g %g \n",                          
            rev.is_inversion, rev.to_acbd, rev.is_fission, rev.same_eds,
            rev.dbreakod[0], rev.dbreakod[1], rev.dbreakod[2], rev.dbreakod[3],
            rev.marker_numbers[0], rev.marker_numbers[1], rev.marker_numbers[2], rev.marker_numbers[3], rev.ld, rev.rd);
}

void
print_rev_brief(FILE* fp, Reversal rev)
{
    fprintf(fp, "%4i %4i %4i %4i, ", (rev.marker_numbers[0]+1)/2, (rev.marker_numbers[1]+1)/2,
            (rev.marker_numbers[2]+1)/2, (rev.marker_numbers[3]+1)/2);    
    if(rev.is_inversion) fprintf(fp, "Inversion. ");
    else fprintf(fp, "Translocation. ");

}

double
f_Ltheat(double T_hot, double T, int Lt)
{
        // grows like (LTHEATFACTOR^Lt) for Lt < Lt_knee,
    
    double x = log(LTHEATFACTOR)/log(T_hot);
    double y1, y2;
        //  printf("T, T_inv, , T_inv^x: %g %g %g \n", T, 1.0/T, pow(T, x));
    y1 = pow(T, (double)(x*(Lt - LTHEATKNEE)));
    y2 = pow(T, (double)(0.1*x*(LTHEATKNEE - Lt)));
        //  y2 = 1.0;

    
    return pow(pow(y1, LTSHARP) + pow(y2, LTSHARP), 1.0/LTSHARP);
 /*    if(Lt < LTHEATKNEE){ */
/*         return pow(T, (double)(x*(Lt - LTHEATKNEE))); */
/*     } */
/*     else{ */
/*         return 1.0; */
/*             // return pow(T, (double)(x*(LTHEATKNEE - Lt))); */
/*     } */
    
        //   return pow((1.0 + pow(alpha, LTHEATFACTOR*LTHEATKNEESHARPNESS*(double)(LTHEATKNEE - Lt))), -1.0/LTHEATKNEESHARPNESS);
}

double f_Ltheat1(double T, double beta, double Lambda, double lambdaT, int Lt)
{     
    double zeta = pow(lambdaT/Lambda, (double)Lt);
  /*   double exponent = (T - pow(T, beta))/pow(T, beta + 1.0); */

/*     return pow(zeta, exponent); */


    return pow(zeta, 1.0/(1.0 + beta*(T - 1.0)))/pow(zeta, 1.0/T);
}


void teststates(Permutation* p)
{
#define ASIZE 200
        int jj;
        double r = 4.0;
        double prob;
        Permutation* perm = copy_perm(p);
        Cycle_decomposition the_cd;
        Permutation* perm_array[ASIZE] = {NULL};
        double perm_count_array[ASIZE] = {0.0};
        double degeneracy_array[ASIZE] = {0.0};
            //    Permutation** p;
        int done = FALSE;
        int ix, iy, lastix = 0;
        int trans_array[ASIZE][ASIZE];
        int nstates = 0;
        double total_rate;
        int M, M0;
        double perm_count_sum = 0.0;
        double degeneracy_sum = 0.0;
        double f, b;
        double x, xmax;
        int kk, nsteps = 0;
        
        printf("state space test \n");
        
        for(kk=0; kk<100; kk++){
        for(jj=0; jj<10000; jj++){
                //  perm = copy_perm(p1[0]); 
            get_CD(perm, perm, &the_cd);
                //   total_rate = perm->n_inv*r + perm->n_trans;
            rand_reverse_r(perm, r, &prob, the_cd.elements);
            total_rate = (double)perm->n_inv*r + (double)perm->n_trans;
            M = perm->n_chrom;
            M0 = M - perm->n_chrom_nonempty;
                //   printf("M, M0: %i %i \n", M, M0);
            for(ix = 0; !done && ix < ASIZE;){
                    //  printf("top of loop over perm_array\n");
                if(perm_array[ix] == NULL){ // nothing in this bin yet, so put perm in it. 
                    perm_array[ix] = copy_perm(perm); // we have to make a copy! since perm will get changed!
                    perm_count_array[ix] = 1.0/total_rate;
                    degeneracy_array[ix] = pow(2.0, (double)(M - M0))*factorial(M)/factorial(M0);
                        //    printf("M, M0, degeneracy: %i %i  %g \n", M, M0, degeneracy_array[ix]);
                        nstates ++;
                    done = TRUE;
                        //   printf("put in new bin, nstates, ix: %i %i \n", nstates, ix);
                }
                else if(are_perms_equivalent(perm_array[ix], perm)){ // this is the place.
                    perm_count_array[ix] += 1.0/total_rate;
                    done = TRUE;
                        //  printf("put in existing bin, nstates, ix, count: %i %i %i\n", nstates, ix, perm_count_array[ix]);
                        // print_perm(stdout, perm); print_perm(stdout, perm_array[ix]);
                }
                else { // this is not the place
                        //  printf("not the right bin, nstates, ix: %i %i \n", nstates, ix);
                    ix++; 
                }     
            }
                //  printf("after loop over perm_array. lastix, ix: %i %i \n", lastix, ix);
            trans_array[lastix][ix]++;
            lastix = ix;
            done = FALSE;
                //  permfree(&perm);
            nsteps++;
        }
    
        xmax = -1.0;
        for(ix=0; ix<nstates; ix++){
            for(iy=0; iy < nstates; iy++){
                f = trans_array[ix][iy]; b = trans_array[iy][ix];
                if(f > 0.0 && b > 0.0){
                x = fabs((f - b)/(f + b));
                if(x > xmax) xmax = x;
                    //  if(f > 0.0 || b > 0.0) printf("%7g \n", (f - b)/(f + b) );
                }
            } //printf("\n");
        }
        printf("nsteps, max imbalance: %i %g \n", nsteps, xmax);
        }

          for(ix=0; ix<nstates; ix++){
            perm_count_sum += perm_count_array[ix];
            degeneracy_sum += degeneracy_array[ix];
        }
        printf("nstates: %i \n", nstates);
        for(ix=0; ix<nstates; ix++){
            printf("%i %g %g \n", perm_array[ix]->n_chrom - perm_array[ix]->n_chrom_nonempty,
                   perm_count_array[ix]/perm_count_sum, degeneracy_array[ix]/degeneracy_sum);
        }
        getchar();
        
    }


int
check_path_revs_ordered(Path* path)
{
    Step* step = path->start;
    int result = TRUE;

    while(step != NULL){
        result = result && check_rev_ordered(step->perm, step->rev);
        step = step->next;
    }
    return result;
}

int rev_just_swaps_chromosome_ends(Permutation* perm, Reversal rev)
{
        // returns TRUE if the argument is a reversal which just swaps 2 chromosome ends
        // resulting in a state which is still member of same equivalence class - i.e.
        // for our purposes, the rev leaves the state unchanged
    int temp; int result;

     if(!rev.to_acbd) {
         assert(!rev.is_inversion);
        temp = rev.marker_numbers[2]; rev.marker_numbers[2] = rev.marker_numbers[3]; rev.marker_numbers[3] = temp;
    }
     if(perm->p[rev.marker_numbers[0]].is_end_marker && perm->p[rev.marker_numbers[3]].is_end_marker) { result = TRUE; }
     else if(perm->p[rev.marker_numbers[1]].is_end_marker && perm->p[rev.marker_numbers[2]].is_end_marker){ result = TRUE; } 
     else result = FALSE;
         //   printf("%8i %8i %8i %8i       %i \n", rev.numbs[0], rev.numbs[1], rev.numbs[2], rev.numbs[3], result);
     return result;
}

int
check_rev_perm_consistent(Permutation* perm, Reversal* rev)
{
        //  Permutation* perm = the_step->perm;
    int a,b,c,d,pa,pb,pc,pd, result = TRUE;
    
    a = rev->marker_numbers[0];
    b = rev->marker_numbers[1];
    c = rev->marker_numbers[2];
    d = rev->marker_numbers[3];

    

    pa = (perm->p[a].is_end_marker)? UNKNOWN: perm->p[a].position;
    pb = (perm->p[b].is_end_marker)? UNKNOWN: perm->p[b].position;
    pc = (perm->p[c].is_end_marker)? UNKNOWN: perm->p[c].position;
    pd = (perm->p[d].is_end_marker)? UNKNOWN: perm->p[d].position;

    if(!(pa == UNKNOWN || pb == UNKNOWN) && (abs(pa-pb) != 1)){
        printf("in check_rev_perm_consistent. a,b positions: %i %i. a,b, numbers: %i %i \n", pa, pb, a, b);
        result = FALSE;
    }
    if(!(pc == UNKNOWN || pd == UNKNOWN) && (abs(pc-pd) != 1)){
        printf("in check_rev_perm_consistent. c,d positions: %i %i. c,d, numbers: %i %i \n", pc, pd, c, d);
        result = FALSE;
    }
        //  if(result == FALSE) getchar();
    return result;
}




// path_etc.c

 
