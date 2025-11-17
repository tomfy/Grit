// perm.c Implementation of permutation handling functions

//#include <iostream>
#include "grit.h"
#include "params.h"
#include "structs.h"
#include "perm.h"
#include "perm_etc.h"
#include "rngs.h"
#include "exts.h"
#include "gr_misc.h"
#include "path_etc.h"

//using namespace std;

#define EPSILON2FACTOR 0.5 // rel probs of steps with delta_c = 1, 0, -1 are 1, epsilon, and epsilon2 = EPSILON2FACTOR*epsilon.


int reverse(Permutation* perm, Reversal* rev)
{
// changes perm to result of applying reversal specified by rev.
    int i, stop, left_chromosome, right_chromosome;
    int left_markers, right_markers, n_inv_lr_old, n_inv_lr_new, delta_n_inv;
    double dlpr = rev->ld + rev->rd;
    marker_end_ptr a, b, c, d, e; // a, b, c, d:  the four ends which will be reconnected
        //  int flip_rh_chrom; // gets set to TRUE in get_abcd if righthand chromosome should be flipped
    double d_ac, d_bd;
    Permutation* perm_in = copy_perm(perm);
    double d_abreak, d_bbreak, d_cbreak, d_dbreak;

    if(!check_rev_ordered(perm, *rev)){
        print_perm(stdout, perm);
        print_rev(stdout, *rev);
        printf("in reverse. check_rev_ordered return FALSE.\n");
        getchar();
    }

#ifdef MYDEBUG
    {
        double ld0, ld1, rd2, rd3;
        ld0 = perm->p[rev->marker_numbers[0]].distance + rev->dbreakod[0]*perm->p[rev->marker_numbers[0]].d_black;
        ld1 = perm->p[rev->marker_numbers[1]].distance - rev->dbreakod[1]*perm->p[rev->marker_numbers[1]].d_black;
        rd2 = perm->p[rev->marker_numbers[2]].distance + rev->dbreakod[2]*perm->p[rev->marker_numbers[2]].d_black;
        rd3 = perm->p[rev->marker_numbers[3]].distance - rev->dbreakod[3]*perm->p[rev->marker_numbers[3]].d_black;

            //    printf("%g %g %g \n", perm->p[rev->marker_numbers[0]].distance, rev->dbreakod[0],  perm->p[rev->marker_numbers[0]].d_black);
            //    printf("%g %g %g \n", perm->p[rev->marker_numbers[1]].distance, rev->dbreakod[1],  perm->p[rev->marker_numbers[1]].d_black);
       
            //     printf("ld0, ld1, rd2,rd3: %g %g %g %g \n", ld0, ld1, rd2, rd3);

            // assert(fabs(ld0-ld1) < 1.0e-12);  assert(fabs(rd2-rd3) < 1.0e-12);
        assert(real_close2(rev->ld, ld0, 1.0e-12, 1.0));
         assert(real_close2(rev->rd, rd2, 1.0e-12, 1.0));
        
      /*   if(!real_close2(rev->ld, ld0, 1.0e-12, 1.0)){ */
/*             printf("rev->ld, rd, ld0, ld1, rd2, rd3: %g %g  %g %g  %g %g \n", rev->ld, rev->rd, ld0, ld1, rd2, rd3);   */
/*             rev->ld = ld0; */
/*             print_perm(stdout, perm); print_rev(stdout, *rev); */
/*             getchar(); */
/*         } */
/*         if(!real_close2(rev->rd, rd2, 1.0e-12, 1.0)){ */
/*             printf("rev->ld, rd, ld0, ld1, rd2, rd3: %g %g  %g %g  %g %g \n", rev->ld, rev->rd, ld0, ld1, rd2, rd3);   */
/*             rev->rd = rd2; */
/*             print_perm(stdout, perm); print_rev(stdout, *rev); getchar(); */
/*         } */
    }
#endif /* MYDEBUG */
    
    a = perm->p + rev->marker_numbers[0];
    b = perm->p + rev->marker_numbers[1];
      c = perm->p + rev->marker_numbers[2];
      d = perm->p + rev->marker_numbers[3]; 
    
// reverse
    if(a->is_end_marker){
        if(!b->is_end_marker && c->is_end_marker && !d->is_end_marker) perm->n_chrom_nonempty--;
        else if(b->is_end_marker && !c->is_end_marker && !d->is_end_marker) perm->n_chrom_nonempty++;
    }
    else{
        if(b->is_end_marker && !c->is_end_marker && d->is_end_marker) perm->n_chrom_nonempty--;
        else if(!b->is_end_marker && c->is_end_marker && d->is_end_marker) perm->n_chrom_nonempty++;
    }
    assert(a->chromosome == b->chromosome);
    assert(c->chromosome == d->chromosome);
    left_chromosome = b->chromosome;
    right_chromosome = c->chromosome;
        // position fields:
        //   printf("ld, rd: %g %g \n", rev->ld, rev->rd);
        //   dlpr = b->distance + c->distance;
       if(rev->ld == UNKNOWN){ // if exact points to break chromosome not specified, just invert so as to preserve distances from a to b, c to d.
            //   printf("ld, rd: %g %g \n", rev.ld, rev.rd);
        assert(rev->rd == UNKNOWN);
        dlpr = b->distance + c->distance;
    }
  
        //   printf("in reverse. rev.ld, rev.rd, dlpr: %g %g %g \n", rev.ld, rev.rd, dlpr);
    if(left_chromosome == right_chromosome){ // inversion 
        for (e = b, stop = b->position, i = c->position; i>= stop; ){
            e->position = i--;
                //    printf("old, new dist: %g  ", e->distance);
            e->distance = dlpr - e->distance; // distances fixed
                //    printf(" %g  \n", e->distance);
            e = e->other;
            e->position = i--;
                //    printf("old, new dist: %g  ", e->distance);
            e->distance = dlpr - e->distance; // distances fixed
                //   printf(" %g  \n", e->distance);
            e = e->black;
        }
    }
    else{ // translocation
        left_markers = perm->n_markers_on_chromosome[left_chromosome]; 
        right_markers = perm->n_markers_on_chromosome[right_chromosome]; 
        n_inv_lr_old = ( (left_markers+1)*left_markers + (right_markers+1)*right_markers )/2;

        for (e = b, stop = b->position, i = c->position; i>= stop; ){
            e->position = i--;
            e->distance = dlpr - e->distance; // distances fixed 
            if(e->chromosome == left_chromosome) {
                e->chromosome = right_chromosome; // swap chromosome numbers in inverted section
                if(!(e->is_end_marker)){
                    perm->n_markers_on_chromosome[left_chromosome]--;
                    perm->n_markers_on_chromosome[right_chromosome]++;
                }
            }
            else if(e->chromosome == right_chromosome) {
                e->chromosome = left_chromosome;
                if(!(e->is_end_marker)){
                    perm->n_markers_on_chromosome[left_chromosome]++;
                    perm->n_markers_on_chromosome[right_chromosome]--;
                }
            }

            e = e->other;
            e->position = i--;
            e->distance = dlpr - e->distance; // distances fixed
            if(e->chromosome == left_chromosome) e->chromosome = right_chromosome;
            else if(e->chromosome == right_chromosome) e->chromosome = left_chromosome;
            e = e->black;
        }
        left_markers = perm->n_markers_on_chromosome[left_chromosome]; 
        right_markers = perm->n_markers_on_chromosome[right_chromosome]; 
        n_inv_lr_new = ( (left_markers+1)*left_markers + (right_markers+1)*right_markers )/2;
        delta_n_inv = n_inv_lr_new - n_inv_lr_old;
        perm->n_inv += delta_n_inv;
        perm->n_trans -= 2*delta_n_inv;

        
      
    } // end of translocation branch
    assert(perm->n_inv == get_n_inv(perm));
 
        // switch the pointers
    a->black = c;
    c->black = a;
    b->black = d;
    d->black = b;

        //   printf("a,b,c,d, d_black: %g %g %g %g \n", a->d_black,  b->d_black,  c->d_black,  d->d_black);
    assert(a->d_black == b->d_black);  assert(c->d_black == d->d_black);
    
    d_abreak = a->d_black*rev->dbreakod[0];
    d_bbreak = b->d_black*rev->dbreakod[1];
    d_cbreak = c->d_black*rev->dbreakod[2];
    d_dbreak = d->d_black*rev->dbreakod[3];
        //   printf("d_abreak, etc: %g %g %g %g \n", d_abreak, d_bbreak,  d_cbreak, d_dbreak); 
    
    d_ac = a->d_black*rev->dbreakod[0] + c->d_black*rev->dbreakod[2];
    d_bd = b->d_black*rev->dbreakod[1] + d->d_black*rev->dbreakod[3];
        //      printf("d_ac, d_bd: %g %g \n", d_ac, d_bd);

    a->d_black = (c->d_black = d_ac);
    b->d_black = (d->d_black = d_bd);
      
        //  assert(check_markend_dists2(perm));
        
    set_Ai_At(perm); // get chromosome lengths, inv/trans areas  
        //  printf("bottom of reverse\n");
    assert(check_perm_selfconsistent(perm));
  /*   if(!check_perm_selfconsistent(perm)) { */
/*         printf("rev dbreakod: %g %g %g %g \n", rev->dbreakod[0], rev->dbreakod[1], rev->dbreakod[2], rev->dbreakod[3]); */
/*         printf("rev rd, ld: %g %g  \n", rev->ld, rev->rd); */
/*         printf("perm in: \n"); print_perm(stdout, perm_in); */
/*         print_perm(stdout, perm); */
/*         getchar(); */
/*     } */
    permfree(&perm_in);
    return TRUE; // success
}  // end of function reverse

Marker_end* get_chromosome_left_end(Marker_end* m)
{
        // given a Marker_end *m,
        // get pointer to Marker_end at left end
        // of chromosome it is on. This will be a chromosome end marker
    
        //  Marker_end* ml;

        // *m should be right end of marker
    assert(m->black->position > m->position);
        // assert(m->other->position < m->position);

    
    for(; !m->is_end_marker; m = m->other->black){
            //  printf("in get_chrom_left_end. m->position %i \n", m->position);
    }
 /*    printf("in get_chrom_left... left end marker: position, is_end_marker: %i %i \n", */
/*            m->position, m->is_end_marker); */
    return m;
}



int reverse_chromosome(Permutation* perm, int the_chromosome, double* rd)
{ // reverses the chromosome specified the_chromosome (1 to n_chrom-1)
        // end markers are reversed along with everything else
        // i.e. "other" edges are broken, not black edges
        // chromosome 0 is never reversed
        // *rd is the 
    int i, stop;
    double dlpr; 
        
    marker_end_ptr the_end, next_end;
    marker_end_ptr a = NULL, b = NULL, c, d, e; // a, b, c, d: the four ends which will be reconnected


        //  printf("top of reverse_chromosome_new\n");
    assert(check_perm_selfconsistent(perm));
        //  if(!check_perm_selfconsistent(perm)) getchar();
    
    the_end = perm->p->black;
    next_end = the_end->other;

    if(the_chromosome <= 0){
        printf("in reverse_chromosome_new, chromosome to reverse: %i \n", the_chromosome);
        return FALSE;
    }
    assert(the_chromosome > 0);
        // find where to cut and reverse; mean # of iterations (2/3)NMARKER
    do {
        if (next_end->chromosome == the_chromosome && the_end->chromosome != the_chromosome){
            a = the_end;
            b = next_end;
        }
        if (next_end->chromosome != the_chromosome && the_end->chromosome == the_chromosome){
            c = the_end;
            d = next_end;
                //  printf("in reverse_chromosome_new, c, d: %p %p \n", c, d);
            break;
        }
            // go to next possible cut site
        the_end = next_end->black;
        next_end = the_end->other;
        if (next_end == NULL) {
            if(the_end->chromosome == the_chromosome){
                c = the_end;
                d = next_end;
                break;
            }
            else{ // should never get here
                fprintf(stdout, "in reverse_chromosome, reached end without finding the chromosome !!!\n");
                process_error();
                return FALSE;
            }
        }
    }
    while(TRUE);

// reverse

    dlpr = b->distance + c->distance;
        //  printf("in reverse_chrom new. *rd, before after: %g ", *rd);
      if(*rd != UNKNOWN){
            //   printf("in reverse_chromosome. rev->ld, rev->rd, before, after: %g %g  ", rev->ld, rev->rd);
        *rd = dlpr - *rd; // now rd is correct for the flipped chromosome
            //  printf("  %g %g \n", rev->ld, rev->rd);
    }
          //    printf(" %g \n", *rd);

    assert(b->chromosome == c->chromosome);
        // position fields:
        //  printf("in reverse_chromosome: dlpr: %g \n", dlpr);
    for (e = b, stop = b->position, i = c->position; i>= stop; ){
        e->position = i--;
            //  printf("old, new dist: %g  ", e->distance);
              e->distance = dlpr - e->distance; // distances fixed
            //  printf("%g  \n", e->distance);
        e = e->black;
        e->position = i--;
            //   printf("old, new dist: %g  ", e->distance);
               e->distance = dlpr - e->distance; // distances fixed
            //  printf("%g  \n", e->distance); 
        e = e->other;
    }
     
        // switch the pointers
    a->other = c;
    c->other = a;
    b->other = d;
    if(d != NULL) d->other = b;
        //  printf("bottom of reverse_chromosome_new\n");
 
        // if(!check_perm_selfconsistent(perm)) getchar();
    assert(check_perm_selfconsistent(perm));
    return TRUE; // success
}  // end of function reverse_chromosome_new

int
check_rev(int n_markers, Reversal* rev)
{
    int mn1 = (rev->marker_numbers[0]+1)/2;
    int mn2 = (rev->marker_numbers[2]+1)/2;
    int checks = TRUE;
    
          if(  ( (mn1 <= n_markers) && (mn1 == (rev->marker_numbers[1]+1)/2) ) ||
               ( (mn2 <= n_markers) && (mn2 == (rev->marker_numbers[3]+1)/2) ) ) {
              printf("rev->marker_numbers: %i %i %i %i \n",
                     rev->marker_numbers[0], rev->marker_numbers[1],
                     rev->marker_numbers[2], rev->marker_numbers[3]);
              checks = FALSE;
                    getchar();
          }
          return checks;
}


void
invert_translocate(Permutation* perm, Reversal* rev, int get_new_dists)
{
        // changes perm to result of applying reversal specified by rev.
        // rev.left, rev.right are marker positions (in range 1 through (n_mark+n_chrom-1))
        //   rather than marker_end positions.
        // rev.toacbd specifies whether to flip the r.h. chromosome involved (ignored if chromosomes are the same)
        // returns TRUE iff rev is an inversion
    
    int flip_it; //, nswaps;
    Reversal revx = *rev;
    Marker_end* A, *B, *C, *D; // the four marker ends which will be reconnected
        //  Marker_end* ends[4], * ordered_ends[4];
        //  double rd = rev->rd;
    double new_A;
    int mn1 = (rev->marker_numbers[0]+1)/2;
    int mn2 = (rev->marker_numbers[2]+1)/2;
        //  int order[4]; // ends[order[i]] is ordered_ends[i]

        //   printf("top of inv...tran..new2.\n");
        //    printf("\n\n top of invert_tran...new. get_new_dists: %i \n", get_new_dists);
          if(  ( (mn1 <= perm->n_mark) && (mn1 == (rev->marker_numbers[1]+1)/2) ) ||
               ( (mn2 <= perm->n_mark) && (mn2 == (rev->marker_numbers[3]+1)/2) ) ) {
              printf("rev->marker_numbers: %i %i %i %i \n",
                     rev->marker_numbers[0], rev->marker_numbers[1],
                     rev->marker_numbers[2], rev->marker_numbers[3]);
              print_perm(stdout, perm);
              getchar();
          }

              //  nswaps = get_4ends(perm, rev->marker_numbers, ends, ordered_ends, order);
    new_A = (perm->p + (rev->marker_numbers[0]))->d_black * (perm->p + (rev->marker_numbers[2]))->d_black;
 
    if(!get_new_dists && !real_close2(rev->A, new_A, 1.0e-12, 1.0)){
        printf("in inv..trans..new. rev->A, new_A: %g %g   %i\n", rev->A, new_A, get_new_dists);
        rev->A = new_A; getchar();
    }
    rev->A = new_A;
        //  flip_it = (nswaps % 2 == 0)? !rev->to_acbd: rev->to_acbd;
        //  flip_it = !rev->to_acbd;
  
    A = perm->p + rev->marker_numbers[0];  B = perm->p + rev->marker_numbers[1];
    C = perm->p + rev->marker_numbers[2];  D = perm->p + rev->marker_numbers[3];
    
         //  A = ordered_ends[0]; B = ordered_ends[1]; C = ordered_ends[2]; D = ordered_ends[3];

    if(C->chromosome != D->chromosome){
        printf("A->chromosome, etc.: %i %i %i %i \n", A->chromosome, B->chromosome, C->chromosome, D->chromosome);
        print_perm(stdout, perm); print_rev(stdout, *rev);
    }
    assert(A->chromosome == B->chromosome); assert(C->chromosome == D->chromosome);
   
        //  is_inversion = (A->chromosome == C->chromosome);
    assert((A->chromosome == C->chromosome) == rev->is_inversion);
    flip_it = (!rev->to_acbd) && (!rev->is_inversion); 
    if(rev->is_inversion) {
        rev->to_acbd = TRUE;
    }
        //   rev->is_inversion = is_inversion;
    rev->is_fission = ((A->is_end_marker && B->is_end_marker) && !(C->is_end_marker || D->is_end_marker))
        || (!(A->is_end_marker || B->is_end_marker) && (C->is_end_marker && D->is_end_marker)); // i.e.
        //  if(rev->is_fission) printf("rev is fission. \n");  
  
    assert(!rev->is_inversion || rev->to_acbd);
// store marker end numbers in numbs array; end_markers all mapped to UNKNOWN
    rev->numbs[0] = A->is_end_marker? UNKNOWN: A->markend_numb;
    rev->numbs[1] = B->is_end_marker? UNKNOWN: B->markend_numb;
    rev->numbs[2] = C->is_end_marker? UNKNOWN: C->markend_numb;
    rev->numbs[3] = D->is_end_marker? UNKNOWN: D->markend_numb;

// only revx (copy of *rev) is modified by reverse_chromosome, then used by reverse
    if(get_new_dists){
        get_rand_break_dists(A, B, C, D, rev->dbreakod, &(rev->ld), &(rev->rd));
    }
    if(flip_it){    
            //  Marker_end* temp_me;
            //  double temp;
        revx = *rev;
          
        reverse_chromosome(perm, D->chromosome, &(revx.rd)); //
            // swap
        revx.dbreakod[2] = rev->dbreakod[3]; revx.dbreakod[3] = rev->dbreakod[2];
        revx.marker_numbers[2] = rev->marker_numbers[3]; revx.marker_numbers[3] = rev->marker_numbers[2];
        
        reverse(perm, &revx);
    }
    else{
        reverse(perm, rev);
    }
       
    check_rev(perm->n_mark, rev);
       
        //   return rev->is_inversion;
} // end of invert_translocate

void
fix_rev1(Permutation* perm, Reversal* rev)
{
        // ignore distances for now at least
    Marker_end* ends[4], * oends[4];
    int i, n_swaps, order[4];
    
    n_swaps = get_4ends(perm, rev->marker_numbers, ends, oends, order);
    for(i=0; i<4; i++){
        rev->marker_numbers[i] = oends[i]->markend_numb;
    }
    if(n_swaps % 2 != 0){
        rev->to_acbd = !rev->to_acbd;
    }

    assert(oends[0]->chromosome == oends[1]->chromosome && oends[2]->chromosome == oends[3]->chromosome);
 
    rev->is_inversion = (oends[0]->chromosome == oends[2]->chromosome);
        //  printf("in fix_rev, check_rev_ordered: %i \n", check_rev_ordered(perm, *rev));
  /*   printf("in fix_rev1. rev->marker_numbers: %i %i %i %i \n", rev->marker_numbers[0], rev->marker_numbers[1], */
/*            rev->marker_numbers[2], rev->marker_numbers[3]); */
           
}
    


int
get_4ends(const Permutation* perm, int numbs[4], Marker_end* ends[4], Marker_end* ordered_ends[4], int order[4])
{   
    int i;
        //   marker_end_ptr  a,b,c,d; // a,b,c,d: the four ends which will be reconnected, a is leftmost in perm, b is a's neighbor (i.e. connected by black edge)
        // d is rightmost in perm and c is d's neighbor
    Marker_end* e01;

        // first check to see whether this is a valid ordered set of ends as is:
    for(i=0; i<4; i++){
        ends[i] = perm->p + numbs[i];
    }
    
    if((ends[0]->black == ends[1]) && (ends[2]->black == ends[3])){ // still valid pairs of neighbors
            //   printf("in get_4ends. ok as is\n");
    }
    else{ // get chrom ends which are actually neighbors of the real markers involved
            //    printf("in get_4ends. not ok as is. fixing\n"); //getchar();
        for(i=0; i<4; i++){
            ends[i] =  perm->p[numbs[i]].is_end_marker? NULL: perm->p + numbs[i];
        }
       
            // find which chrom ends are there in perm, but don't put in ends in L to R order
        e01 = fix_me_pair(perm, ends, NULL);
        (void)fix_me_pair(perm,  ends+2, e01);
    }
    if(!real_close2(ends[0]->d_black, ends[1]->d_black, 1.0e-12, 1.0) || !real_close2(ends[2]->d_black, ends[3]->d_black, 1.0e-12, 1.0)){
        printf("ends[0]->d_black, etc: %g %g  %g %g \n", ends[0]->d_black, ends[1]->d_black, ends[2]->d_black, ends[3]->d_black); 
        printf("ends[0], [1], ends[0]->black, ends[1]->black: %p %p %p %p \n", ends[0], ends[1], ends[0]->black, ends[1]->black);
        printf("ends[2], [3], ends[2]->black, ends[3]->black: %p %p %p %p \n", ends[2], ends[3], ends[2]->black, ends[3]->black);
        printf("ends[0,1]->markend_numb %i  %i  \n", ends[0]->markend_numb, ends[1]->markend_numb);
        printf("ends[0,1]->black->markend_numb,  %i  %i  \n", ends[0]->black->markend_numb, ends[1]->black->markend_numb);
        printf("ends[2,3]->markend_numb %i  %i  \n", ends[2]->markend_numb, ends[3]->markend_numb);
        printf("ends[2,3]->black->markend_numb,  %i  %i  \n", ends[2]->black->markend_numb, ends[3]->black->markend_numb);
            //   print_perm(stdout, perm);
        getchar();
    }
     assert((ends[0]->d_black == ends[1]->d_black) && (ends[2]->d_black == ends[3]->d_black));
    return order_4ends(ends, ordered_ends, order); 
}  // get_4_ends

Marker_end* 
fix_me_pair(const Permutation* perm, Marker_end** e, Marker_end* prev_e_empty)
{
        // n is array of 2 marker_end numbers of marker_ends separated by black edge
        // i.e. which are adjacent in perm
        // pos is array of their positions, UNKNOWN if a marker_end is a chrom_end
        // puts marker_ends in order s.t. e[0]->position +1 = e[1]->position
        // finds correct chrom_end to use if one chrom end and one real marker, by taking
        // marker_end adjacent (in perm) to the real one as the chrom end
        // if empty chromosome, (i.e. both e[0], e[1] are chrom ends) returns ptr to one end of the empty chrom; else returns NULL
#define FIXMEPAIRNEW TRUE
    Marker_end* e_empty = NULL;
        //  int nne = 2*(perm->n_mark + perm->n_chrom_nonempty), n = 2*(perm->n_mark + perm->n_chrom);
    int nme; 
    if(e[0] != NULL){
        if(e[1] != NULL){
        }
        else{
            e[1] = e[0]->black; //
            assert(e[1]->black == e[0]);
        }
    }
    else{ // e[0] is NULL
        if(e[1] != NULL){
            e[0] = e[1]->black; //
            assert(e[1] == e[0]->black);
        }
        else{ // both e[0], e[1] are NULL - both chrom ends - markerless chromosome
                //  e_empty = e[0];

            if(FIXMEPAIRNEW){ // choose a chrom end marker at random
                do{
                    nme = (int)(2*perm->n_chrom*drand(&rng_long)); 
                    e[0] = (nme == 0)? perm->p: perm->p + 2*perm->n_mark + nme;
                    assert(e[0]->is_end_marker);
                    e[1] = e[0]->black;
                    
                }while(!(e[1]->is_end_marker && (e[0] != prev_e_empty && e[1] != prev_e_empty)));
                e_empty = e[0];
            }
            else{
                
                e[0] = perm->p+(2*perm->n_mark);
                do{
                    e[0] = e[0]++;
                    assert(e[0]->is_end_marker);
                    e[1] = e[0]->black; 
                }while(!(e[1]->is_end_marker && (e[0] != prev_e_empty && e[1] != prev_e_empty)));
                e_empty = e[0];
            }
        }
    }
        //   printf("in fix_me_pair. prev_e_empty, e[0], e[1]: %i %i %i \n", prev_e_empty, e[0], e[1]);
    /*   if(prev_e_empty != NULL) printf("in fix_me_pair. prev_e_empty->position, e[0]->position, e[1]->position: %i %i %i \n", */
/*            prev_e_empty->position, e[0]->position, e[1]->position); */
/*       if(e_empty != NULL) printf("in fix_me_pair. e_empty: %i  (should be NULL)\n", e_empty); */
/*       printf("prev_e_empty, e_empty: %i %i; e[0], e[1]: %i %i \n", prev_e_empty, e_empty, e[0], e[1]); */
    return e_empty;
} // end of fix_me_pair

int
order_4ends(Marker_end* ends[4], Marker_end* ordered_ends[4], int order[4])
{
        // returns number of swaps (0, 1, or 2)
    Marker_end* temp;
    int nswaps = 0;
    int i;
    
        // ends is array of 4 Marker_ends*, at end ordered_ends
        // has same pointers, but in L to R order
        //   printf("ends[0]->position, ends[1]->position: %i %i \n", ends[0]->position, ends[1]->position);
    assert(ends[0]->position/2 == ends[1]->position/2);
    if(ends[0]->position < ends[1]->position){
        ordered_ends[0] = ends[0]; ordered_ends[1] = ends[1];
        order[0] = 0; order[1] = 1;
    }
    else if(ends[0]->position > ends[1]->position){
        ordered_ends[0] = ends[1]; ordered_ends[1] = ends[0];
        order[0] = 1; order[1] = 0;
        nswaps++;
    }
    else printf("in order_4ends. two end positions are same. \n");
    
    assert(ends[2]->position/2 == ends[3]->position/2);
    if(ends[2]->position < ends[3]->position){
        ordered_ends[2] = ends[2]; ordered_ends[3] = ends[3];
        order[2] = 2; order[3] = 3;
    }
    else if(ends[2]->position > ends[3]->position){
        ordered_ends[2] = ends[3]; ordered_ends[3] = ends[2];
        order[2] = 3; order[3] = 2;
        nswaps++;
    }
    else printf("in order_4ends. two end positions are same. \n");

  if(ordered_ends[0]->position > ordered_ends[2]->position){ // swap first two and last two
      temp = ordered_ends[0]; ordered_ends[0] = ordered_ends[2]; ordered_ends[2] = temp;
      temp = ordered_ends[1]; ordered_ends[1] = ordered_ends[3]; ordered_ends[3] = temp;
      int_swap(order, order+2); int_swap(order+1, order+3);
  }
  for(i=0; i<4; i++){
      if(!(ordered_ends[i] == ends[order[i]])){
          printf("oend[i]->distance, ends[order[i]]->distance: %g %g \n",
                 ordered_ends[i]->distance, ends[order[i]]->distance); getchar();
      }
  }
  return nswaps;
} // end of order_4ends

Reversal rand_reverse_r(Permutation* perm, double r, double* prob, Cycle_element* elements)
{
        // choose one of the (N+M choose 2) pairs of breakpoints
        // do this inversion or translocation, with, in the case
        // of a translation, a 1/2 probability of flipping one
        // of the chromosomes first
        // r is ratio of prob of generating a given inv to prob of generating a given translocation.
        // if r <= 2 accept inversions with prob r/2
        // if r > 2 accept translocations with prob 2/r
    int n_edges = 0;
    int e1, e2;
    Reversal the_rev;
    int is_inversion, done = FALSE;
        //  Marker_end* a, * b, * c, * d;
    int flipped;
    
    n_edges = perm->n_mark + perm->n_chrom;

    if(elements == NULL){
        printf("in rand_reverse_r. Called with elements == NULL\n"); process_error();
    }
    do{
        e1 = (int)((n_edges)*drand(&rng_long));
        e2 = (int)((n_edges-1)*drand(&rng_long));
            //  printf("in rand_reverse_r. n_edges, e1, e2: %i %i %i \n", n_edges, e1, e2);
   
        if (e2 >= e1) e2++;
            //  flipped = (int)(2*drand(&rng_long));
            
        is_inversion = (elements[e1].chromosome == elements[e2].chromosome);
        flipped = (is_inversion)? FALSE: (int)(2*drand(&rng_long));
            
            //    if(r<2.0) done = !is_inversion || (drand(&rng_long) >= 0.5*r);
            //  if(r<2.0) done = ((!is_inversion) || (drand(&rng_long) >= 0.5*r));
        if(r<2.0) done = ((!is_inversion) || (drand(&rng_long) < 0.5*r));
        else done = is_inversion || drand(&rng_long) < 2.0/r;
            
    }while(!done ); //
      
    the_rev = make_rev(n_edges, elements + e1, elements + e2, flipped, perm);
        //  assert((the_rev.left == e1+1 && the_rev.right == e2) || (the_rev.left == e2+1 && the_rev.right == e1));
            
    assert(get_n_inv(perm) == get_n_inv1(perm) && perm->n_inv == get_n_inv(perm));
    *prob = is_inversion? r: 1.0;
    *prob *= 1.0/((double)perm->n_trans + r*(double)perm->n_inv);
    invert_translocate(perm, &the_rev, FALSE);
    assert(the_rev.is_inversion == is_inversion);
    return the_rev;
}  // end of function rand_reverse_r

inline marker_end_ptr get_cycle_start(marker_end_ptr the_end, int* cycle_found)
{
        // move L to R through marker ends until find one which is not yet on a cycle, return ptr to that one.
        // assumes the_end is at the left end of a black edge 
   
     while((the_end != NULL) && (cycle_found[the_end->markend_numb] == TRUE)) {         
        assert((the_end->black->position - the_end->position) == 1); // the_end should be at the left end of a black edge
        the_end = the_end->black->other;
    }
    return the_end;
} // end of function get_cycle_start

Reversal make_signflip_rev(Permutation* perm, int cut1, int cut2)
{
  // cut1, cut2 are the positions of black edges to be cut
  // cut1, cut2 should both be in range 0 to n_edges-1
        // distance ld is chosen uar from (d1, d2) where d1, d2 are distances of markers bounding the cut edge
        // and similarly for rd

    int i, temp, done = FALSE;
    Reversal the_result = r_init;
    Marker_end* end, * next_end, *a = NULL, *b = NULL, *c = NULL, *d = NULL;
    int n_edges = perm->n_mark + perm->n_chrom;
        // double min_sep;

    the_result.to_acbd = TRUE;
    
    if ((cut1 < 0) || (cut1 >= n_edges) || (cut2 < 0) || (cut2 >= n_edges)) {
		fprintf(stdout, "In make_signflip_rev, cut 1 or cut 2 out of range. \n");
      printf("nedges, cut1, cut2: %i %i %i \n", n_edges, cut1, cut2);
      process_error();
		return the_result;
	}

   if(cut1 > cut2){ // swap so cut1 < cut2
       temp = cut1; cut1 = cut2; cut2 = temp;
   }
   assert(cut1+1 ==  cut2);
  
   end = perm->p;
   next_end = end->black; // end, next_end are adjacent in perm; connected by black edge
   do{
       assert(end->position/2 == next_end->position/2);
       if(end->position/2 == cut1){
           the_result.marker_numbers[0] = end->markend_numb;
           the_result.marker_numbers[1] = next_end->markend_numb;
           a = end; b = next_end;
               //  assert(a->d_black == b->d_black);
           if(!(a->d_black == b->d_black)){
                 printf("a,b d_black %g %g ; a.b numb: %i %i \n", a->d_black, b->d_black, a->markend_numb, b->markend_numb);
               print_perm(stdout, perm);  getchar();
           }
       }
       
       if(end->position/2 == cut2){ 
           the_result.marker_numbers[2] = end->markend_numb;
           the_result.marker_numbers[3] = next_end->markend_numb;
           c = end; d = next_end;
               //  assert(c->d_black == d->d_black);
           if(!(c->d_black == d->d_black)){
                printf("c,d d_black %g %g ; c,d numb: %i %i \n", c->d_black, d->d_black, c->markend_numb, d->markend_numb);
                print_perm(stdout, perm); getchar();
           }
           done = TRUE;
       }
       else{
           end = next_end->other;
           next_end = end->black;
       }
   }while(!done);

   
       // set dbreakod so as to leave ab and cd separation distance unchanged by flip (as appropriate for sign flips!)
   get_signflip_break_dists(a, b, c, d, the_result.dbreakod, &(the_result.ld), &(the_result.rd));

   the_result.is_inversion = a->chromosome == d->chromosome;

   for(i=0; i<4; i++){
           the_result.numbs[i] = perm->p[the_result.marker_numbers[i]].is_end_marker? UNKNOWN: the_result.marker_numbers[i];
   }

       if(!check_rev_ordered(perm, the_result))
           printf("in make_signflip_rev\n");
  assert( check_rev_ordered(perm, the_result));
   
	return the_result;
}  // end of function make_signflip_rev

void
get_signflip_break_dists(Marker_end* a, Marker_end* b, Marker_end* c, Marker_end* d, double dbreakod[4], double* ld, double* rd)
{
    double min_sep = (a->d_black < c->d_black)? a->d_black: c->d_black;
  
    dbreakod[1] = 0.5*min_sep/b->d_black;
    dbreakod[2] = 0.5*min_sep/c->d_black;
    dbreakod[0] = 1.0 - dbreakod[1];
    dbreakod[3] = 1.0 - dbreakod[2];
    *ld = a->distance*dbreakod[1] + b->distance*dbreakod[0];
    *rd = c->distance*dbreakod[3] + d->distance*dbreakod[2];
}

void
get_rand_break_dists(Marker_end* a, Marker_end* b, Marker_end* c, Marker_end* d, double dbreakod[4], double* ld, double* rd)
{
    double rand;
    
    rand = drand(&rng_long);
    dbreakod[0] = rand;
    dbreakod[1] = 1.0 - rand;
    *ld = a->distance*dbreakod[1] + b->distance*dbreakod[0];
    
    rand = drand(&rng_long); 
    dbreakod[2] = rand;
    dbreakod[3] = 1.0 - rand;
    *rd = c->distance*dbreakod[3] + d->distance*dbreakod[2];
}

Reversal make_rev1(const Permutation* perm, Marker_end* a, Marker_end* b, Marker_end* c, Marker_end* d,
                   int flipped_chromosome)
{
    Reversal the_result = r_init;

    the_result.marker_numbers[0] = a->markend_numb;
    the_result.marker_numbers[1] = b->markend_numb;
    the_result.marker_numbers[2] = c->markend_numb;
    the_result.marker_numbers[3] = d->markend_numb;

  get_rand_break_dists(a, b, c, d, the_result.dbreakod, &(the_result.ld), &(the_result.rd));
    
    the_result.is_inversion = (a->chromosome == c->chromosome);
    the_result.to_acbd = (flipped_chromosome)? FALSE: TRUE;

    the_result.A = a->d_black * c->d_black;
    assert(real_close2(a->d_black, fabs(a->distance - b->distance), 1.0e-12, 1.0));
    assert(real_close2(c->d_black, fabs(c->distance - d->distance), 1.0e-12, 1.0));

    return the_result;

} // end make_rev1

void simple_rand_reverse(Permutation* perm, double r)
{
  // choose one of the (N+M choose 2) pairs of breakpoints
        // do this inversion or translocation, with, in the case
        // of a translation, a 1/2 probability of flipping one
        // of the chromosomes first
        // r is ratio of prob of generating a given inv to prob of generating a given translocation.
        // if r <= 2 accept inversions with prob r/2
        // if r > 2 accept translocations with prob 2/r
        // like rand_reverse_r, except doesn't require Cycle_element* argument
        // doesn't return reversal or calculate prob
        // primarily for use in randomizing marker order initially
    int n_edges = 0;
    int e1, e2;
    Reversal the_rev;
    int is_inversion, done = FALSE;
    Marker_end* a = NULL, * b = NULL, * c = NULL, * d = NULL, *ml, *mr;
    int flipped;
    int etmp, i_edge;
    
    n_edges = perm->n_mark + perm->n_chrom;
    do{
        e1 = (int)((n_edges)*drand(&rng_long));
        e2 = (int)((n_edges-1)*drand(&rng_long));
            //  printf("in rand_reverse_r. n_edges, e1, e2: %i %i %i \n", n_edges, e1, e2);
   
        if (e2 >= e1) e2++;
        if(e1 > e2){ etmp = e1; e1 = e2; e2 = etmp; } //swap so that e2 > e1
            //  flipped = (int)(2*drand(&rng_long));
            // run along permutation, find a,b,c,d
        for(i_edge = 0, ml = perm->p, mr = ml->black; ml != NULL; ml = mr->other, i_edge++){
            mr = ml->black;
            if(i_edge == e1){ a = ml; b = mr; }
            if(i_edge == e2){ c = ml; d = mr; break; }
        }
            
        is_inversion = (a->chromosome == c->chromosome);
        flipped = (is_inversion)? FALSE: (int)(2*drand(&rng_long));
            
            //    if(r<2.0) done = !is_inversion || (drand(&rng_long) >= 0.5*r);
            //  if(r<2.0) done = ((!is_inversion) || (drand(&rng_long) >= 0.5*r));
        if(r<2.0) done = ((!is_inversion) || (drand(&rng_long) < 0.5*r));
        else done = is_inversion || drand(&rng_long) < 2.0/r;
            
    }while(!done ); //
      
    the_rev = make_rev1(perm, a, b, c, d, flipped);
        //  assert((the_rev.left == e1+1 && the_rev.right == e2) || (the_rev.left == e2+1 && the_rev.right == e1));
            
    assert(get_n_inv(perm) == get_n_inv1(perm) && perm->n_inv == get_n_inv(perm));
        //  *prob = is_inversion? r: 1.0;
        // *prob *= 1.0/((double)perm->n_trans + r*(double)perm->n_inv);
    invert_translocate(perm, &the_rev, FALSE);
    assert(the_rev.is_inversion == is_inversion);
        // return the_rev;
}  // end of function simple_rand_reverse


Reversal make_rev(int n_edges, Cycle_element* edge1, Cycle_element* edge2, int flip_chromosome, const Permutation* perm)
{
  // edge1, edge2 are the black edges to be cut
        // n_edges = (n_mark + n_chrom)
        // flip_chromosome specifies whether to flip the right chromosome (ignored if boths cuts are on one chromosome)
        // distance ld is chosen uar from (d1, d2) where d1, d2 are distances of markers bounding the cut edge
        // and similarly for rd

        // if a subpath and it's replacement subpath are not identical (i.e. although in the same equivalence class, they
        // differ by chromosome end relabelling, chromosome flipping and swapping)
        // then when we set new_subpath->end->rev = old_subpath->end->rev in update path, the marker_numbers (and positions)
        // will not in general be correct to give the desired next perm after the subpath; in particular chromosome ends
        // may be wrong. So don't use marker_numbers of chromosome ends - get them from 

    Reversal the_result = r_init;
       
    int cut1 = edge1->position;
    int cut2 = edge2->position;
    double d1, d2;
        //  double A = 1.0;
    double the_rand_number;
    Marker_end* a, * b, * c, * d;
    
	if ((cut1 < 0) || (cut1 >= n_edges) || (cut2 < 0) || (cut2 >= n_edges)) {
		fprintf(stdout, "In make_rev, cut 1 or cut 2 out of range. \n");
      printf("nedges, cut1, cut2: %i %i %i \n", n_edges, cut1, cut2);
      process_error();
		return the_result;
	}
   
	if (cut1 < cut2) {
           // set marker_numbers; 0,1 at ends of left edge; 2,3 right,
           // and 0 to left of 1, 2 to left of 3 .. so they are in order 0-3 L to R
       if(edge1->direction == RIGHT){
           assert(edge1->d2 > edge1->d1);
           assert(perm->p[edge1->n1].position + 1 == perm->p[edge1->n2].position);
           assert(perm->p[edge1->n1].position/2 == cut1);
           the_result.marker_numbers[0] = edge1->n1; the_result.marker_numbers[1] = edge1->n2;
       }
       else{
           assert(edge1->d2 < edge1->d1);
            assert(perm->p[edge1->n1].position - 1 == perm->p[edge1->n2].position);
           assert(perm->p[edge1->n1].position/2 == cut1);
           the_result.marker_numbers[0] = edge1->n2; the_result.marker_numbers[1] = edge1->n1;
       }
       if(edge2->direction == RIGHT){
            assert(perm->p[edge2->n1].position + 1 == perm->p[edge2->n2].position);
           assert(perm->p[edge2->n1].position/2 == cut2);
           assert(edge2->d2 > edge2->d1);
           the_result.marker_numbers[2] = edge2->n1; the_result.marker_numbers[3] = edge2->n2;
       }
       else{
            assert(perm->p[edge2->n1].position - 1 == perm->p[edge2->n2].position);
           assert(perm->p[edge2->n1].position/2 == cut2); 
           assert(edge2->d2 < edge2->d1);
           the_result.marker_numbers[2] = edge2->n2; the_result.marker_numbers[3] = edge2->n1;
       }
      
    /*    the_result.left = cut1 + 1; */
/*        the_result.right = cut2; */
       
       d1 = edge1->d1; d2 = edge1->d2;
       if(d1 > d2) double_swap(&d1, &d2);
       the_rand_number = drand(&rng_long);
       the_result.dbreakod[0] = the_rand_number;
       the_result.dbreakod[1] = (1.0 - the_rand_number);
       
       the_result.ld = d1 + the_rand_number*(d2 - d1);
       if(DBDIST && (d1 < 0.0 || d2 < 0.0)) {printf("in make_rev. d1, d2: %g %g \n", d1, d2);getchar();}
         
              assert(real_close2(edge1->d12, fabs(d2 - d1), 1.0e-11, 1.0));
       
       d1 = edge2->d1; d2 = edge2->d2;
        if(d1 > d2) double_swap(&d1, &d2);
        if(DBDIST && (d1 < 0.0 || d2 < 0.0)) {printf("in make_rev. d1, d2: %g %g \n", d1, d2);getchar();}
       the_rand_number = drand(&rng_long);
       the_result.dbreakod[2] = the_rand_number;
       the_result.dbreakod[3] = (1.0 - the_rand_number);
       
       the_result.rd = d1 + the_rand_number*(d2 - d1);
             assert(real_close2(edge2->d12, fabs(d2 - d1), 1.0e-11, 1.0));
   }
   else if (cut1 > cut2) {
       if(edge2->direction == RIGHT){
           the_result.marker_numbers[0] = edge2->n1; the_result.marker_numbers[1] = edge2->n2;
           assert(edge2->d2 > edge2->d1);
           assert(perm->p[edge2->n1].position + 1 == perm->p[edge2->n2].position);
           assert(perm->p[edge2->n1].position/2 == cut2);
       }
       else{
           assert(edge2->d2 < edge2->d1);
           assert(perm->p[edge2->n1].position - 1 == perm->p[edge2->n2].position);
           assert(perm->p[edge2->n1].position/2 == cut2);
           the_result.marker_numbers[0] = edge2->n2; the_result.marker_numbers[1] = edge2->n1;
       }
       if(edge1->direction == RIGHT){
           assert(edge1->d2 > edge1->d1);
           assert(perm->p[edge1->n1].position + 1 == perm->p[edge1->n2].position);
           assert(perm->p[edge1->n1].position/2 == cut1);
           the_result.marker_numbers[2] = edge1->n1; the_result.marker_numbers[3] = edge1->n2;
       }
       else{
           assert(edge1->d2 < edge1->d1);
           assert(perm->p[edge1->n1].position - 1 == perm->p[edge1->n2].position);
           assert(perm->p[edge1->n1].position/2 == cut1);
           the_result.marker_numbers[2] = edge1->n2; the_result.marker_numbers[3] = edge1->n1;
       }
       
       d1 = edge1->d1; d2 = edge1->d2;
        if(d1 > d2) double_swap(&d1, &d2);
       if(DBDIST && (d1 < 0.0 || d2 < 0.0)) {printf("in make_rev. d1, d2: %g %g \n", d1, d2);getchar();}
       the_rand_number = drand(&rng_long);
       the_result.dbreakod[2] = the_rand_number;
       the_result.dbreakod[3] = (1.0 - the_rand_number);
       the_result.rd = d1 + the_rand_number*(d2 - d1);
         
       assert(real_close2(edge1->d12, fabs(d2 - d1), 1.0e-11, 1.0)); 
       d1 = edge2->d1; d2 = edge2->d2;
       if(d1 > d2) double_swap(&d1, &d2);
       if(DBDIST && (d1 < 0.0 || d2 < 0.0)) {printf("in make_rev. d1, d2: %g %g \n", d1, d2);getchar();}
       the_rand_number = drand(&rng_long);
       the_result.dbreakod[0] = the_rand_number;
       the_result.dbreakod[1] = (1.0 - the_rand_number);
       the_result.ld = d1 + the_rand_number*(d2 - d1);
         
       assert(real_close2(edge2->d12, fabs(d2 - d1), 1.0e-11, 1.0));
              
	}
   else {
       fprintf(stdout, "In make_rev, cut1 == cut2.\n");
       process_error();
	}
   the_result.A = edge1->d12*edge2->d12;
   the_result.to_acbd = !flip_chromosome;
      
   if(! check_rev_ordered(perm, the_result)){
       print_perm(stdout, perm); print_rev(stdout, the_result); getchar(); 
   }

   a = perm->p + the_result.marker_numbers[0];
   b = perm->p + the_result.marker_numbers[1];
   c = perm->p + the_result.marker_numbers[2];
   d = perm->p + the_result.marker_numbers[3];
   assert(a->chromosome == b->chromosome && c->chromosome == d->chromosome);
 
   the_result.is_inversion = (a->chromosome == c->chromosome);
   the_result.n_chroms_involved = the_result.is_inversion? 1: 2;
   the_result.chroms_involved[0] = a->chromosome;
   if(the_result.n_chroms_involved == 2) the_result.chroms_involved[1] = c->chromosome;
   if(the_result.marker_numbers[0] == the_result.marker_numbers[2] || the_result.marker_numbers[0] == the_result.marker_numbers[3]){
       print_rev(stdout, the_result);
   }
   return the_result;
}  // end of function make_rev

Permutation* copy_perm(const Permutation* src)
{
    int i, imax = 2*(src->n_mark + src->n_chrom);
  Marker_end*  sptr = src->p,* tptr;
  Marker_end* t, * s;
  Permutation* targ = multipermalloc(src->n_chrom, src->n_mark);
  
  tptr = targ->p;
  *targ = *src; // copy fields of Permutation
  targ->p = tptr;
  t = tptr; s = sptr;
 
  // fix the pointers
  for (i = 0; i < imax; i++, t++, s++){
          
      *t = *s; // set each marker end equal to corresponding one in src, but then need to
          // use markend_numb (= array index) of objects pointed to in src to fix black, other ptrs
        
      t->black = tptr + s->black->markend_numb;  // fix the pointers.
  
      if(t->other != NULL) t->other = tptr + s->other->markend_numb;  // fix the pointers.
    
  }
  assert(check_perm_selfconsistent(targ));

  return targ;
} // end of function copy_perm





int get_CD(const Permutation* perm1, const Permutation* perm2, Cycle_decomposition* the_cd)
{
        // new version, nov 15 02
    int the_cycle_number = 0;
    marker_end_ptr cycle_start = perm1->p; // ptr to marker end to start a cycle
    int n_chromosome_ends = 0;
    
    { // initialize the_cd
        int i; int* c;
        the_cd->perm1 = perm1; the_cd->perm2 = perm2;
        the_cd->n_mark = perm1->n_mark;
        the_cd->n_chrom = perm1->n_chrom;
        the_cd->n_empty_chrom_in_target = 0;
        for(i=0; i<5; i++) {
            the_cd->n_dsw_inv_ddm1[i] =  the_cd->n_dsw_trans_ddm1[i] = 0;
            the_cd->n_ddm1_nsw2before[i] = 0;
            the_cd->n_ddm1_ncomdectozero[i] = 0;
            the_cd->n_ddm1_ncomdecrease[i] = 0;
            the_cd->n_ddm1_ncomnochange[i] = 0;
        }
       
            // set all cycle_found[i] to FALSE
        for(i=0, c=the_cd->cycle_found; i<2*(perm1->n_mark + perm1->n_chrom); i++, c++){ *c = FALSE; }
        the_cd->n_end_cycles = 0;
        the_cd->n_edss = 0;
        for(i=0; i<3; i++){
            the_cd->n_dci[i] = 0;
            the_cd->n_dd[i] = 0; // not actually needed
            the_cd->n_dcg[i] = 0; // not actually needed
            the_cd->n_dd_inv[i] = 0;
            the_cd->n_dd_trans[i] = 0;
        }
        for(i=0; i<(perm1->n_mark + perm1->n_chrom); i++){
            the_cd->elements[i].position = i; // elements are in position order
            the_cd->elements[i].next = NULL;
        }
    }  // end of initialization

//    printf("in get_CD, before loop\n");
    while ((cycle_start = get_cycle_start(cycle_start, the_cd->cycle_found)) != NULL){
            //  printf("in get_CD, in loop. cycle_state: %p \n", cycle_start);
        n_chromosome_ends += run_around_cycle(perm1, perm2, the_cd, cycle_start);
        the_cycle_number++;
    } // end of while loop over cycles
//printf("in get_CD, after loop\n");

    assert(n_chromosome_ends == 2*perm1->n_chrom); 
    the_cd->n_cycles = the_cycle_number;
    the_cd->n_int_cycles = the_cd->n_cycles - the_cd->n_end_cycles;
     
  
        // get n_dcim1 from total minus the other 4 possibilities
        //  if(ddgezero_possible(the_cd) == FALSE) {printf("ddgezero_possible returned FALSE!\n"); getchar();}
        //  assert(ddgezero(the_cd));  
    get_end_int_edss(the_cd);
             
    the_cd->n_dci[2] = get_n_dci1(the_cd, perm2);
    the_cd->n_dcg[0] = get_n_dcgm1(the_cd, perm2);
  
    the_cd->n_dd[0] = the_cd->n_dcg[0] + the_cd->n_dci[2];
    the_cd->n_dd_inv[0] = the_cd->n_dcg_inv[0] + the_cd->n_dci_inv[2];
    the_cd->n_dd_trans[0] = the_cd->n_dcg_trans[0] + the_cd->n_dci_trans[2];

    the_cd->c_g = get_c_g(the_cd);
    the_cd->d = the_cd->n_mark - the_cd->n_chrom - the_cd->n_int_cycles + the_cd->c_g; // 
        //   printf("d: %i \n", the_cd->d);
        // used to get all the n_dd[i] here  
  
    assert(perm2->n_chrom - perm2->n_chrom_nonempty == the_cd->n_empty_chrom_in_target);
    assert(the_cd->n_int_cycles == (the_cd->n_edss + the_cd->n_empty_chrom_in_target - 2*(the_cd->n_chrom)));
      //   check_edss(the_cd);
   
    return the_cycle_number;
} // end of function get_CD


int
run_around_cycle(const Permutation* perm1, const Permutation* perm2,  Cycle_decomposition* the_cd,
                 marker_end_ptr cycle_start)
{ // this version runs around the cycles of the breakpoint graph of perm1 relative to perm2
  // i.e. its black edges are the adjacencies of perm1, its gray edges are the adjacencies of perm2
        // returns number of chromosome ends on the cycle

        // new version, oct 28 02
    marker_end_ptr the_next_end, the_end = cycle_start; // initialize the_end to starting point (marker_end) of this cycle
    Cycle_element* ce;
    Cycle_element* the_elem_ptr, * eds_start_ptr;
    int n_ends_on_cycle = 0; // counts the number of chromosome end markers on the cycle
    int eds_number_in_cycle = 0; // number of the end_delimited_segment
    int start_eds_number = the_cd->n_edss;
      
    eds_start_ptr = the_cd->elements + the_end->position/2; // i.e. pointer to edge position/2 (0 to N+M-1)
    the_cd->edss[the_cd->n_edss].start_elem_ptr = eds_start_ptr;
    the_cd->edss[the_cd->n_edss].is_end_cycle = TRUE; // set to FALSE later if no ends on cycle (and therefore only 1 eds)
    the_cd->edss[the_cd->n_edss].n_black = 0;
    the_cd->edss[the_cd->n_edss].n_gray = 0;
       
    do {
        the_next_end = the_end->black;
        
        the_cd->cycle_found[the_end->markend_numb]
            = the_cd->cycle_found[the_next_end->markend_numb] = TRUE;  // record that marker_end's cycle is now known
        
        assert(the_end->position/2 == the_next_end->position/2); // the 2 positions should be an even number and its successor
        the_elem_ptr = the_cd->elements + the_end->position/2;
        assert(the_end->position/2 == the_elem_ptr->position);
        
        if (the_next_end->position == the_end->position + 1) {  // black edge is crossed rightward
            the_elem_ptr->direction =  RIGHT;
            the_elem_ptr->at_left_end = the_end->is_end_marker;
            the_elem_ptr->at_right_end = the_next_end->is_end_marker;
        }
        else { // (the_next_end->position == the_end->position - 1);  black edge is crossed leftward
            the_elem_ptr->direction = LEFT;
            the_elem_ptr->at_left_end = the_next_end->is_end_marker;
            the_elem_ptr->at_right_end = the_end->is_end_marker;
        }

        the_elem_ptr->position = the_end->position/2; // position of the black edge
        the_elem_ptr->n1 = the_end->markend_numb;
        the_elem_ptr->n2 = the_next_end->markend_numb;
        the_elem_ptr->d1 = the_end->distance; // distance from L end of genomes
        the_elem_ptr->d2 = the_next_end->distance;
        the_elem_ptr->d12 = the_end->d_black;
      /*   if(!real_close2(the_end->d_black, the_next_end->d_black, 1.0e-12, 1.0)){ */
/*             printf("in run_around_cycle. the_end->d_black, the_next_end->d_black: %g %g \n", the_end->d_black, the_next_end->d_black); */
/*             printf("in run_around_cycle. the_end->black->d_black, the_next_end->black->d_black: %g %g \n", the_end->black->d_black, the_next_end->black->d_black); */
/*             print_perm(stdout, perm1); */
/*         } */
        assert(the_end->d_black == the_next_end->d_black);
        the_elem_ptr->chromosome = the_end->chromosome;
        assert(the_end->chromosome == the_next_end->chromosome);
        
        the_elem_ptr->eds = the_cd->n_edss;
        
            // follow one black edge and one gray edge
            // i.e. go to next marker end on perm->eds_in_cycle1, and then to the marker end which is next to that on perm2<<
        the_end = perm1->p + ((perm2->p[the_next_end->markend_numb].black)->markend_numb);
            // now the_end, the_next_end joined by gray edge (i.e. adjacent in perm2)
        the_elem_ptr->next = the_cd->elements + the_end->position/2;

        the_cd->edss[the_cd->n_edss].n_black++;
        if(the_next_end->is_end_marker) {
            the_cd->edss[the_cd->n_edss].last_black = TRUE;  
            n_ends_on_cycle++;
            eds_number_in_cycle++; // new eds
            the_cd->n_edss++;
            the_cd->edss[the_cd->n_edss].is_end_cycle = TRUE;
            the_cd->edss[the_cd->n_edss].n_black = 0; 
            eds_start_ptr = the_cd->elements + the_end->position/2; // store ptr to this edge which is first of a new eds
            if(the_end->is_end_marker) {
                    // there is an eds with gray on both ends (indeed just a single gray edge, empty chromosome in target genome (perm2))
                    // skip this and go to the next cycle, which begins with the_end
                the_cd->edss[the_cd->n_edss].first_black = TRUE;
                the_cd->edss[the_cd->n_edss].is_end_cycle = TRUE;
                the_cd->edss[the_cd->n_edss].n_gray = 0;
                n_ends_on_cycle++;       
                the_cd->n_empty_chrom_in_target++; // successive marker_ends separated by gray edge are end markers. -> empty chrom in perm2
            }
            else{
                the_cd->edss[the_cd->n_edss].first_black = FALSE;
                the_cd->edss[the_cd->n_edss].is_end_cycle = TRUE;
                the_cd->edss[the_cd->n_edss].n_gray = 1;
            }
        }
        else { // the_next_end is not an end marker
            the_cd->edss[the_cd->n_edss].n_gray++;
            if (the_end->is_end_marker) {
                the_cd->edss[the_cd->n_edss].last_black = FALSE;  
                n_ends_on_cycle++;
                eds_number_in_cycle++; // new eds
                the_cd->n_edss++;
                the_cd->edss[the_cd->n_edss].first_black = TRUE;
                the_cd->edss[the_cd->n_edss].is_end_cycle = TRUE;
                the_cd->edss[the_cd->n_edss].n_black = 0;
                the_cd->edss[the_cd->n_edss].n_gray = 0; 
                eds_start_ptr = the_cd->elements + the_end->position/2; // store ptr to this edge which is first of a new eds
            }
        }
        the_cd->edss[the_cd->n_edss].start_elem_ptr = eds_start_ptr;
    } while (the_end != cycle_start);
       
    if(eds_number_in_cycle > 0){ // end cycle
        the_cd->n_end_cycles++;
    }
    else{ // int cycle
        the_cd->edss[the_cd->n_edss].first_black = TRUE; // so int cycles have first_black == TRUE
        the_cd->edss[the_cd->n_edss].last_black = FALSE; // so int cycles have last_black == FALSE
        the_cd->edss[the_cd->n_edss].is_end_cycle = FALSE;  //
        the_cd->n_edss++;
    }
    
    if(eds_number_in_cycle > 0){ // end cycle
            // set cycle to start at beginning of an eds
        the_cd->edss[start_eds_number].start_elem_ptr = eds_start_ptr;
        the_cd->edss[start_eds_number].first_black = the_cd->edss[the_cd->n_edss].first_black;
        the_cd->edss[start_eds_number].n_black += the_cd->edss[the_cd->n_edss].n_black;
        the_cd->edss[start_eds_number].n_gray += the_cd->edss[the_cd->n_edss].n_gray;
     
            // set eds of this (last) eds on cycle to eds number of first eds on cycle
        ce = eds_start_ptr;
        for(ce = eds_start_ptr; ce->eds == the_cd->n_edss; ce = ce->next){
            ce->eds = start_eds_number;
        }
    }
    return n_ends_on_cycle;
} // end of function run_around_cycle


int
get_c_g(Cycle_decomposition* the_cd)
{
    Cycle* the_eds;
    int i, cg1 = 0, cg2 = 0;
    for(i=0; i<the_cd->n_edss; i++){
        the_eds = the_cd->edss + i;
        if(the_eds->is_end_cycle){
            if(the_eds->first_black && the_eds->last_black) cg1++;
            else if(!the_eds->first_black && !the_eds->last_black) cg2++;
        }
    }
    cg2 += the_cd->n_empty_chrom_in_target;
    assert(cg1 == cg2);
    return cg1;
}



int
get_n00_disteds2(Cycle_decomposition* the_cd)
{
        // looks at pairs of distinct edss
        // of type bb+bb, bb+bg, gg+gg, etc. i.e. with number of black eds ends != 2;
        // these all give delta c_i = delta c_g = 0
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j, nb1;
    int n_00 = 0;
    int n00_deds2_inv = 0, n00_deds2_trans = 0;
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle){
            nb1 = eds1->first_black + eds1->last_black;
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle && ((nb1 + eds2->first_black + eds2->last_black) != 2)){ 
                    edge1 = eds1->start_elem_ptr;
                    do{
                        edge2 = eds2->start_elem_ptr;
                        do{
                            if(edge1->chromosome == edge2->chromosome){ n_00++; n00_deds2_inv++;}
                            else { n_00 += 2; n00_deds2_trans += 2;}
                            edge2 = edge2->next;
                        }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                        edge1 = edge1->next;
                    }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                }
            }
        }
    }
    the_cd->n00_deds2_inv = n00_deds2_inv;
    the_cd->n00_deds2_trans = n00_deds2_trans;
    return n_00;
} // end of get_n00_disteds2

Reversal
get_spec_dcgm1_rev(Cycle_decomposition* the_cd, int n, int type)
{
        // returns the number of delta c_g = -1 inversion/translocations
    Cycle* eds1, * eds2;
    Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcgm1 = 0,  n_dcgm1_inv = 0, n_dcgm1_trans = 0 ;
    int* n_dcgm1_type;

    if(type == EITHER){ n_dcgm1_type = &n_dcgm1; }
    else if(type == INVERSION){ n_dcgm1_type = &n_dcgm1_inv; }
    else if(type == TRANSLOCATION){ n_dcgm1_type = &n_dcgm1_trans; }
    else{ n_dcgm1_type = NULL; printf("in get_spec_00seds_rev. unknown type. \n");}
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle && (eds1->first_black == eds1->last_black)){
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle  && (eds2->first_black == eds2->last_black) && (eds1->first_black != eds2->first_black)){
                    edge1 = eds1->start_elem_ptr;
                    do{
                        edge2 = eds2->start_elem_ptr;
                        do{
                            if(edge1->chromosome == edge2->chromosome){ //inversion
                                n_dcgm1++; n_dcgm1_inv++;
                                if(*n_dcgm1_type == n){ // inversion
                                    return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                }
                            }
                            else { // translocation
                                n_dcgm1 += 2; n_dcgm1_trans += 2;
                                if(*n_dcgm1_type == n+1){ // unflipped translocation
                                    return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                }
                                else if (*n_dcgm1_type == n){ // flipped translocation
                                    return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                                }
                            }
                            edge2 = edge2->next;
                        }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                        edge1 = edge1->next;
                    }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                }
            }
        }
    }
    printf("in get_spec_dcgm1_rev. Reached end without finding specified reversal.\n");
} // end of get_spec_dcgm1_rev


int
get_n_dcim1(Cycle_decomposition* the_cd)
{
        // returns the number of delta c_i = -1 inversion/translocations
        // uses edss, not cycles
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcim1 = 0;
    int n_dcim1_inv = 0, n_dcim1_trans = 0;
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
            if((!eds1->is_end_cycle) || (!eds2->is_end_cycle)){
                edge1 = eds1->start_elem_ptr;
                do{
                    edge2 = eds2->start_elem_ptr;
                    do{
                        if(edge1->chromosome == edge2->chromosome){ n_dcim1++; n_dcim1_inv++;}
                        else { n_dcim1 += 2; n_dcim1_trans += 2;}
                        edge2 = edge2->next;
                    }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                    edge1 = edge1->next;
                }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
            }
        }
    }
    the_cd->n_dci_inv[0] = n_dcim1_inv;
    the_cd->n_dci_trans[0] = n_dcim1_trans;
    return n_dcim1;
} // end of n_dcim1


        

Reversal
get_spec_dcim1_rev(Cycle_decomposition* the_cd, int n, int type)
{
        // returns the n_th reversal (inversion/translocations)
        // with delta c_i = -1
        // uses edss, not cycles
        // type is INV, TRANS or EITHER
    Cycle* eds1, * eds2;
    Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcim1 = 0, n_dcim1_inv = 0, n_dcim1_trans = 0;
    int* n_dcim1_type;

    if(type == EITHER){ n_dcim1_type = &n_dcim1; }
    else if(type == INVERSION){ n_dcim1_type = &n_dcim1_inv; }
    else if(type == TRANSLOCATION){ n_dcim1_type = &n_dcim1_trans; }
    else{ n_dcim1_type = NULL; printf("in get_spec_dcim1_rev. unknown type. \n");}
        
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
            if((!eds1->is_end_cycle) || (!eds2->is_end_cycle)){
                edge1 = eds1->start_elem_ptr;
                do{
                    edge2 = eds2->start_elem_ptr;
                    do{
                        if(edge1->chromosome == edge2->chromosome){ // inversion
                            n_dcim1++; n_dcim1_inv++;
                            if(n == *n_dcim1_type) return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                        }
                        else { // translocation
                            n_dcim1 += 2; n_dcim1_trans += 2;
                            if(*n_dcim1_type == n+1) return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                            else if(*n_dcim1_type == n) return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                        }
                        edge2 = edge2->next;
                    }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                    edge1 = edge1->next;
                }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
            }
        }
    }
    printf("In get_spec_dcim1_rev. Reached end without finding specified rev.\n");
} // end of get_spec_dcim1_rev




int count_between_chrom_breaks(const Cycle_element* edge1, const Cycle_element* edge2, const Permutation* p2)
{
    
    int nnnn = 0;
    int chrom_a, chrom_b, chrom_c, chrom_d;
    Marker_end* a, * b, * c, * d;    
        //
    a = p2->p + edge1->n1;
    b = p2->p + edge1->n2;
    c = p2->p + edge2->n1;
    d = p2->p + edge2->n2;
                            
    chrom_a = a->chromosome;
    chrom_b = b->chromosome;
    chrom_c = c->chromosome;
    chrom_d = d->chromosome;  
    if(!(a->is_end_marker || b->is_end_marker) && chrom_a != chrom_b) nnnn++;
    if(!(c->is_end_marker || d->is_end_marker) && chrom_c != chrom_d) nnnn++;
    return nnnn;
}

int n_syn(const Permutation* p1, const Permutation * p2)
{
        // runs along p1, counting changes in the chromosome on p2
    Marker_end* e, * next_e;
    int switches = 0;

    for(e = p1->p; e != NULL;  e = next_e->other){
        next_e = e->black;
        if(!p1->p[e->markend_numb].is_end_marker && !p1->p[next_e->markend_numb].is_end_marker){
            switches += (p2->p[e->markend_numb].chromosome == p2->p[next_e->markend_numb].chromosome)? 0: 1;
        }
        else switches++; // at least one marker end
    }

    return switches;
}
    

int delta_n_switches(const Cycle_element* edge1, const Cycle_element* edge2, const Permutation* p2, int* flipped_result,
                     int* nswbefore, int* nswafter, int* nswafterflipped)
{
        // imagine that the markers are colored according to which chromosome of p2 they are on, and there is an additional color
        // which is given to all chromosome ends
        // How many times in going from one end of genome to other does the color change; call this number n_switches?
        // n_switches is a measure of how mixed up the chromosomes are.
        // We want to prefer steps which reduce n_switches. This subroutine calculates the change in n_switches from the reversal involving
        // edge1 and edge2. The return value is delta_n_switches for inversions or unflipped translocation. flipped result points to n_switches
        // for a flipped translocation
        // int nnnn = 0;
    int chrom_a, chrom_b, chrom_c, chrom_d;
    Marker_end* a, * b, * c, * d;
    int switches_before = 0, switches_after = 0, switches_after_flipped = 0;
        //  int pa, pb, pc, pd;

    if(edge1->direction == LEFT){
        a = p2->p + edge1->n1;
        b = p2->p + edge1->n2;
    }
    else{
        a = p2->p + edge1->n2;
        b = p2->p + edge1->n1;
    }
    if(edge2->direction == LEFT){
        c = p2->p + edge2->n1;
        d = p2->p + edge2->n2;
    }
    else{
        c = p2->p + edge2->n2;
        d = p2->p + edge2->n1;
    }

    if(CHROMENDMATCHING == SELF){ 
    chrom_a = (a->is_end_marker)? UNKNOWN: a->chromosome;
    chrom_b = (b->is_end_marker)? UNKNOWN: b->chromosome;
    chrom_c = (c->is_end_marker)? UNKNOWN: c->chromosome;
    chrom_d = (d->is_end_marker)? UNKNOWN: d->chromosome;
    }
    else if(CHROMENDMATCHING == NONE){  // chrom ends match nothing
        chrom_a = (a->is_end_marker)? -1: a->chromosome;
        chrom_b = (b->is_end_marker)? -2: b->chromosome;
        chrom_c = (c->is_end_marker)? -3: c->chromosome;
        chrom_d = (d->is_end_marker)? -4: d->chromosome;
    }
    else printf("In delta_n_switches. chrom end matching option %i not implemented \n", CHROMENDMATCHING);
    switches_before += (chrom_a == chrom_b)? 0:1;
    switches_before += (chrom_c == chrom_d)? 0:1;
    
    switches_after += (chrom_a == chrom_c)? 0:1;
    switches_after += (chrom_b == chrom_d)? 0:1;

    switches_after_flipped += (chrom_a == chrom_d)? 0:1;
    switches_after_flipped += (chrom_b == chrom_c)? 0:1;
 /*    printf("in delta_n_switches. chroms a,b,c,d: %i %i %i %i    nsw, bef, aft, aftflip: %i %i %i \n", */
/*            chrom_a, chrom_b, chrom_c, chrom_d, switches_before, switches_after, switches_after_flipped); */
    *nswbefore = switches_before;
    *nswafter = switches_after;
    *nswafterflipped = switches_after_flipped;
    
    *flipped_result = switches_after_flipped - switches_before;
    return switches_after - switches_before;;
}

int delta_n_switches1(const int* marker_numbs, const Permutation* p1, const Permutation* p2, int* flipped_result)
{
        // like delta_n_switches, except that just give array of marker_end numbers ("marker_numbs") instead
        // of getting them from edge1->n1, etc.
        // imagine that the markers are colored according to which chromosome of p2 they are on, and there is an additional color
        // which is given to all chromosome ends
        // How many times in going from one end of genome to other does the color change; call this number n_switches?
        // n_switches is a measure of how mixed up the chromosomes are.
        // We want to prefer steps which reduce n_switches. This subroutine calculates the change in n_switches from the reversal involving
        // edge1 and edge2. The return value is delta_n_switches for inversions or unflipped translocation. flipped result points to n_switches
        // for a flipped translocation
        //  int nnnn = 0;
    int chrom_a, chrom_b, chrom_c, chrom_d;
    Marker_end* a = NULL, * b = NULL , * c = NULL, * d = NULL;
    int switches_before = 0, switches_after = 0, switches_after_flipped = 0;
    int pa, pb, pc, pd;
    
    pa = p1->p[marker_numbs[0]].position;
    pb = p1->p[marker_numbs[1]].position;
    pc = p1->p[marker_numbs[2]].position;
    pd = p1->p[marker_numbs[3]].position;
    
    if(pb == pa+1){
        a = p2->p + marker_numbs[0];
        b = p2->p + marker_numbs[1];
    }
    else if(pa == pb+1){
        a = p2->p + marker_numbs[1];
        b = p2->p + marker_numbs[0];
    }
    else {
        printf("a,b not in neighboring positions. positions of a,b: %i %i .marker_number: %i %i  \n",
               pa, pb, marker_numbs[0], marker_numbs[1]);
        print_perm(stdout, p1);
    }
    if(pd == pc+1){    
        c = p2->p + marker_numbs[2];
        d = p2->p + marker_numbs[3];
    }
    else if(pc == pd+1){    
        c = p2->p + marker_numbs[3];
        d = p2->p + marker_numbs[2];
    }
    else  {
        printf("c,d not in neighboring positions: positions of c,d: %i %i   .marker_number: %i %i \n",
               pc, pd, marker_numbs[2], marker_numbs[3]);
        print_perm(stdout, p1);
    }

  
    if(CHROMENDMATCHING == SELF){  // chrom ends match other chrom ends
        chrom_a = (a->is_end_marker)? UNKNOWN: a->chromosome;
        chrom_b = (b->is_end_marker)? UNKNOWN: b->chromosome;
        chrom_c = (c->is_end_marker)? UNKNOWN: c->chromosome;
        chrom_d = (d->is_end_marker)? UNKNOWN: d->chromosome;
    }
    else if(CHROMENDMATCHING == NONE){  // chrom ends match nothing
        chrom_a = (a->is_end_marker)? -1: a->chromosome;
        chrom_b = (b->is_end_marker)? -2: b->chromosome;
        chrom_c = (c->is_end_marker)? -3: c->chromosome;
        chrom_d = (d->is_end_marker)? -4: d->chromosome;
    }
        //   printf("in ..sw1. chroms, a,b,c,d: %i %i %i %i \n", chrom_a, chrom_b, chrom_c, chrom_d);
    
    switches_before += (chrom_a == chrom_b)? 0:1;
    switches_before += (chrom_c == chrom_d)? 0:1;
    
    switches_after += (chrom_a == chrom_c)? 0:1;
    switches_after += (chrom_b == chrom_d)? 0:1;

    switches_after_flipped += (chrom_a == chrom_d)? 0:1;
    switches_after_flipped += (chrom_b == chrom_c)? 0:1;
    
    *flipped_result = switches_after_flipped - switches_before;
    return switches_after - switches_before;;
}

int delta_n_switches2(const int* marker_numbers, const Permutation* p1, const Permutation* p2, int* flipped_result)
{
        // like delta_n_switches, except that just give array of marker_end numbers ("marker_numbs") instead
        // of getting them from edge1->n1, etc.
        // imagine that the markers are colored according to which chromosome of p2 they are on, and there is an additional color
        // which is given to all chromosome ends
        // How many times in going from one end of genome to other does the color change; call this number n_switches?
        // n_switches is a measure of how mixed up the chromosomes are.
        // We want to prefer steps which reduce n_switches. This subroutine calculates the change in n_switches from the reversal involving
        // edge1 and edge2. The return value is delta_n_switches for inversions or unflipped translocation. flipped result points to n_switches
        // for a flipped translocation
        //  
    int chrom_a, chrom_b, chrom_c, chrom_d;
        //  Marker_end* a, * b, * c, * d;
    int switches_before = 0, switches_after = 0, switches_after_flipped = 0;
        //int pa, pb, pc, pd;
    int i;
    int numbs[4];

          printf("in ...sw2. ");
          printf("marker_numbers\n");
    for(i=0; i<4; i++){
                  printf(" %i ", marker_numbers[i]);
        numbs[i] = (p1->p[marker_numbers[i]].is_end_marker)? UNKNOWN: marker_numbers[i];
    }printf("\n");

         printf("positions #s\n");
    for(i=0; i<4; i++){
        printf(" %i ", p1->p[marker_numbers[i]].position);
    }printf("\n"); 
    printf("chromosome #s\n");
    for(i=0; i<4; i++){
        printf(" %i ", p2->p[marker_numbers[i]].chromosome);
    }printf("\n");
    
    get_chroms_for_switches(numbs[0], numbs[1], TRUE, &chrom_a, &chrom_b, p1, p2);
    get_chroms_for_switches(numbs[2], numbs[3], FALSE, &chrom_c, &chrom_d, p1, p2);
           printf("in ..sw2. chroms, a,b,c,d: %i %i %i %i \n", chrom_a, chrom_b, chrom_c, chrom_d);
    switches_before += (chrom_a == chrom_b)? 0:1;
    switches_before += (chrom_c == chrom_d)? 0:1;
    
    switches_after += (chrom_a == chrom_c)? 0:1;
    switches_after += (chrom_b == chrom_d)? 0:1;

    switches_after_flipped += (chrom_a == chrom_d)? 0:1;
    switches_after_flipped += (chrom_b == chrom_c)? 0:1;
    
    *flipped_result = switches_after_flipped - switches_before;
    return switches_after - switches_before;;
} // end delta_n_switches2

int n_switches(const int* marker_numbers, const Permutation* perm, const Permutation* targ_perm)
{
        // get the number of syntenic segment switches between a-b and c-d, (0, 1, or 2)
    int i, switches = 0;
    int numbs[4];
    int chrom_a, chrom_b, chrom_c, chrom_d;
    
    for(i=0; i<4; i++){ numbs[i] = (perm->p[marker_numbers[i]].is_end_marker)? UNKNOWN: marker_numbers[i]; }
    
    get_chroms_for_switches(numbs[0], numbs[1], TRUE, &chrom_a, &chrom_b, perm, targ_perm);
    get_chroms_for_switches(numbs[2], numbs[3], FALSE, &chrom_c, &chrom_d, perm, targ_perm);
         
    switches += (chrom_a == chrom_b)? 0:1;
    switches += (chrom_c == chrom_d)? 0:1;

    return switches;
}

int n_switches1(const int* marker_numbers, const Permutation* perm, const Permutation* targ_perm)
{
        // get the number of syntenic segment switches at the two edges determined by marker_numbers
        // assumes that either 02 share an edge (and 13), or 03 and 12.    
    int i, matches = 0;
    int numbs[4], chroms[4], positions[4];
        // int chrom_a, chrom_b, chrom_c, chrom_d;
    int real_02, real_13, real_03, real_12, real_01, real_23;
    
    for(i=0; i<4; i++){
       
        numbs[i] = (perm->p[marker_numbers[i]].is_end_marker)? UNKNOWN: marker_numbers[i];
        chroms[i] = (numbs[i] == UNKNOWN)? -(i+1): targ_perm->p[numbs[i]].chromosome;
        positions[i] = (numbs[i] == UNKNOWN)? -10*(i+1): perm->p[marker_numbers[i]].position;
        printf(" %i ", numbs[i]);
    }printf("\n");

    real_02 = (numbs[0] != UNKNOWN && numbs[2] != UNKNOWN);
    real_13 = (numbs[1] != UNKNOWN && numbs[3] != UNKNOWN);
    real_03 = (numbs[0] != UNKNOWN && numbs[3] != UNKNOWN);
    real_12 = (numbs[1] != UNKNOWN && numbs[2] != UNKNOWN);
    real_01 = (numbs[0] != UNKNOWN && numbs[1] != UNKNOWN);
    real_23 = (numbs[2] != UNKNOWN && numbs[3] != UNKNOWN);

    if( (real_01  && (positions[0]/2 == positions[1]/2)) || (real_23  && (positions[2]/2 == positions[3]/2)) ){ //ab cd
        printf("n_switches1, ab cd branch\n");
        if(real_01) matches += (chroms[0] == chroms[1])? 1: 0;
        if(real_23) matches += (chroms[2] == chroms[3])? 1: 0;
    }
    else if( (real_02  && (positions[0]/2 == positions[2]/2)) || (real_13  && (positions[1]/2 == positions[3]/2)) ){ //ac bd
         printf("n_switches1, ac bd branch\n");
        if(real_02) matches += (chroms[0] == chroms[2])? 1: 0;
        if(real_13) matches += (chroms[1] == chroms[3])? 1: 0;
    }
    else if( (real_03  && (positions[0]/2 == positions[3]/2)) || (real_12  && (positions[1]/2 == positions[2]/2)) ){ //ad bc
         printf("n_switches1, ad cb branch\n");
        if(real_03) matches += (chroms[0] == chroms[3])? 1: 0;
        if(real_12) matches += (chroms[1] == chroms[2])? 1: 0;
    }
    else { // no matches
    }
    
    return 2 - matches;
}

void
get_chroms_for_switches(int n0, int n1, int ab, int* chrom_a, int* chrom_b, const Permutation* p1, const Permutation* p2)
{
    int pos0, pos1;
    int offset = (ab)? 0: -2; //

        // look at n0, and n1, which are either marker_end numbers representing real (neighboring)
        // marker_ends, or UNKNOWN representing a chromosome end
        // *chrom_a, *chrom_b, are the chromosome numbers (in genome p2) of the left and right of n1, n2;
        // but chromosome ends are assigned chromosome numbers which are all distinct from each other
        // offset is 0 for a,b, 2 for c,d. 
    if(n0 != UNKNOWN){
        pos0 = p1->p[n0].position;
        if(n1 != UNKNOWN){ // usual case - both a,b are real
            pos1 = p1->p[n1].position;
            if(pos1 == pos0+1){ //no swap
                assert(pos0 % 2 == 0);
                *chrom_a = p2->p[n0].chromosome;
                *chrom_b = p2->p[n1].chromosome;
            }
            else if(pos0 == pos1+1){ //swap
                assert(pos0 % 2 == 1);
                *chrom_a = p2->p[n1].chromosome;
                *chrom_b = p2->p[n0].chromosome;
            }
            else {
                printf("In get_chroms_for_switches. positions or a,b are not neighboring\n");
                printf("a,b, numbers, positions: %i %i    %i %i\n", n0, n1, pos0, pos1);
                print_perm(stdout, p1); print_perm(stdout, p2);
            }
        }
        else{ // n1 - chromosome end
            if(pos0 %2 == 0){ // no swap
                *chrom_a = p2->p[n0].chromosome;
                *chrom_b = -2 + offset;
            }
            else{
                *chrom_a = -1 + offset;
                *chrom_b = p2->p[n0].chromosome;
            }
        }
    }
    else{ // n0 == UNKNOWN
        if(n1 != UNKNOWN){
            pos1 =  p1->p[n1].position;
            if(pos1 %2 == 0){ // swap
                *chrom_a = p2->p[n1].chromosome;
                *chrom_b = -2 + offset;
            }
            else{
                *chrom_a = -1 + offset;
                 *chrom_b = p2->p[n1].chromosome;
            }
        }
        else{ // a,b both unknown
            *chrom_a = -1 + offset;
            *chrom_b = -2 + offset;
        }
    }
} 

void
get_end_int_edss(Cycle_decomposition* the_cd)
{
        // store pointers to all the end edss in array the_cd->end_edss;
    int i, n_end_edss = 0, n_int_edss = 0;
    Cycle* eds;
  
    for(i=0, eds = the_cd->edss; i<the_cd->n_edss; i++, eds++){
        if(eds->is_end_cycle){
            the_cd->end_edss[n_end_edss] = eds;
            n_end_edss++;
        }
        else{
           the_cd->int_edss[n_int_edss] = eds;
            n_int_edss++;
        } 
    }
    the_cd->n_end_edss = n_end_edss;
    assert(the_cd->n_end_edss + the_cd->n_empty_chrom_in_target == 2*the_cd->n_chrom);
    assert(the_cd->n_int_cycles == n_int_edss);
} // end get_end_edss


int
get_n_dcg1_old(Cycle_decomposition* the_cd)
{
        // returns the number of delta c_g = 1 inversion/translocations
    Cycle* eds1, * eds2;
    Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcg1 = 0;
    int n = 0;
        
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle && (eds1->first_black != eds1->last_black)){
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle  && (eds2->first_black != eds2->last_black)){ 
                    if(eds1->first_black == eds2->first_black){ // bg + bg or gb + gb, same (opp) directions required (unflipped) for delta c_d = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        do{
                            edge2 = eds2->start_elem_ptr;
                            do{
                                if(edge1->chromosome == edge2->chromosome){
                                    if(edge1->direction == edge2->direction) n_dcg1++;
                                    else n++;
                                }
                                else { n_dcg1 += 1; n++; } // exactly 1 of flipped, unflipped will be delta c_g = 1
                                edge2 = edge2->next;
                            }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                            edge1 = edge1->next;
                        }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                    }
                    else { // bg + gb, opposite (same) directions required (unflipped) for delta c_d = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        do{
                            edge2 = eds2->start_elem_ptr;
                            do{
                                if(edge1->chromosome == edge2->chromosome){
                                    if(edge1->direction != edge2->direction) n_dcg1++;
                                    else n++;
                                }
                                else { n_dcg1 += 1; n++; } // exactly 1 of flipped, unflipped will be delta c_g = 1
                                edge2 = edge2->next;
                            }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                            edge1 = edge1->next;
                        }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                    }
                }
            }
        }
    }
    the_cd->n00_deds1 = n;
    return n_dcg1;
} // end of get_n_dcg1_old

int
get_n_dcg1(Cycle_decomposition* the_cd)
{
    // returns the number of delta c_g = 1 inversion/translocations
    // these all arise from pairs of distinct edss of type bg+bg, bg+gb, etc.
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcg1 = 0;
    int n_dcg1_inv = 0, n_dcg1_trans = 0;
    int n00_deds1_inv = 0, n00_deds1_trans = 0;
    int n = 0;
    int ii, jj;
        
    for(i=0, eds1 = the_cd->end_edss[0]; i<the_cd->n_end_edss; i++, eds1 = the_cd->end_edss[i]){
        assert(eds1->is_end_cycle == TRUE);
        if(eds1->first_black != eds1->last_black){
            for(j = i+1, eds2 = the_cd->end_edss[j]; j < the_cd->n_end_edss; j++, eds2 = the_cd->end_edss[j]){
                if(eds2->first_black != eds2->last_black){
                    assert(eds2->is_end_cycle == TRUE);
                    if(eds1->first_black == eds2->first_black){ // bg + bg or gb + gb, same (opp) directions required (unflipped) for delta c_d = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        for(ii=0; ii<eds1->n_black; ii++, edge1=edge1->next){
                            edge2 = eds2->start_elem_ptr;
                            for(jj=0; jj<eds2->n_black; jj++, edge2=edge2->next){
                                if(edge1->chromosome == edge2->chromosome){ // inversions
                                    if(edge1->direction == edge2->direction){ n_dcg1++; n_dcg1_inv++; }
                                    else {n++; n00_deds1_inv++;}
                                }
                                else { n_dcg1 += 1; n_dcg1_trans++;  n++; n00_deds1_trans++;} // exactly 1 of flipped, unflipped will be delta c_g = 1
                            }
                        }
                    }
                    else { // bg + gb, opposite (same) directions required (unflipped) for delta c_d = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        for(ii=0; ii<eds1->n_black; ii++, edge1=edge1->next){
                            edge2 = eds2->start_elem_ptr;
                            for(jj=0; jj<eds2->n_black; jj++, edge2=edge2->next){
                                if(edge1->chromosome == edge2->chromosome){
                                    if(edge1->direction != edge2->direction) {n_dcg1++; n_dcg1_inv++; }
                                    else {n++; n00_deds1_inv++;}
                                }
                                else { n_dcg1 += 1; n_dcg1_trans++; n++;  n00_deds1_trans++;} // exactly 1 of flipped, unflipped will be delta c_g = 1
                            }
                        }
                    }
                }
            }
        }
    }
    the_cd->n00_deds1 = n;
    the_cd->n_dcg_inv[2] = n_dcg1_inv;
    the_cd->n_dcg_trans[2] = n_dcg1_trans;
    the_cd->n00_deds1_inv = n00_deds1_inv;
    the_cd->n00_deds1_trans = n00_deds1_trans; 
    return n_dcg1;
} // end of get_n_dcg1



Reversal
get_spec_dcg1_rev(Cycle_decomposition* the_cd, int n, int type)
{
    // returns the n_th Reversal (inversion/translocation) with delta c_g = 1
    // type is one of EITHER, INVERSION, TRANSLOCATION
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcg1 = 0, n_dcg1_inv = 0, n_dcg1_trans = 0;
    int* n_dcg1_type;

    if(type == EITHER){ n_dcg1_type = &n_dcg1; }
    else if(type == INVERSION){ n_dcg1_type = &n_dcg1_inv; }
    else if(type == TRANSLOCATION){ n_dcg1_type = &n_dcg1_trans; }
    else{ n_dcg1_type = NULL; printf("in get_spec_dcg1_rev. unknown type. \n");}
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle && (eds1->first_black != eds1->last_black)){
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle  && (eds2->first_black != eds2->last_black)){ 
                    if(eds1->first_black == eds2->first_black){ // bg + bg or gb + gb, same (opp) directions required (unflipped) for delta c_d = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        do{
                            edge2 = eds2->start_elem_ptr;
                            do{
                                if(edge1->chromosome == edge2->chromosome){ //inversion
                                    if(edge1->direction == edge2->direction){
                                        n_dcg1++; n_dcg1_inv++;
                                        if(*n_dcg1_type == n){
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                    }
                                }
                                else {
                                    n_dcg1 += 1; n_dcg1_trans++;// exactly 1 of flipped, unflipped will be delta c_g = 1
                                    if(*n_dcg1_type == n){
                                        if(edge1->direction == edge2->direction){
                                           return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                        else{
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                                        }
                                    }
                                }
                                edge2 = edge2->next;
                            }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                            edge1 = edge1->next;
                        }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                    }
                    else { // bg + gb, opposite (same) directions required (unflipped) for delta c_d = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        do{
                            edge2 = eds2->start_elem_ptr;
                            do{
                               if(edge1->chromosome == edge2->chromosome){ //inversion
                                    if(edge1->direction != edge2->direction){
                                        n_dcg1++; n_dcg1_inv++;
                                        if(*n_dcg1_type == n){
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                    }
                                }
                                else {
                                    n_dcg1 += 1; n_dcg1_trans++; // exactly 1 of flipped, unflipped will be delta c_g = 1
                                    if(*n_dcg1_type == n){
                                        if(edge1->direction != edge2->direction){
                                           return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                        else{
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                                        }
                                    }
                                }
                                edge2 = edge2->next;
                            }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                            edge1 = edge1->next;
                        }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                    }
                }
            }
        }
    }
    printf("In get_spec_dcg1_rev. Reached end without finding specified reversal.\n");
} // end of get_spec_dcg1_rev


int
get_genomes_from_file(FILE* fp, Permutation** genomes, int Use_distances)
{
    int zzz = 43;
   
        // if Use_distances is TRUE and
        // if first genome in file has distances uses them, and
        // spaces markers in second genome equally along genome, with
        // total length same as in first genome
        // if no distances, just use L_g = N+m for both, equally spaced (i.e. spacing is 1.0)
        // return value is TRUE if at least one genome has distance information
    int Total_N_markers = -1, N_genomes = -1, N_chromosomes;
    int N_markers[2][MAX_N_CHROMOSOMES] = {{0}};
    Marker perm_array[2][MAX_N_CHROMOSOMES][MAX_N_MARKERS_PER_CHROMOSOME];
    int i,j,k, i_genome, i_chromosome;
    Permutation* the_perm;
    Marker g[2][NMARKERMAX];
    int i_marker, n_marker[2] = {0,0};
    int distances_in_file = FALSE, distances_this_genome; // TRUE/FALSE, indicates whether file has marker distances
    double dist;
    char tempstr[128], temp2[128];
    double chromosome_length[2][MAX_N_CHROMOSOMES];
    double cumulative_length_of_chromosomes = 0.0;
    double L_g, marker_spacing;
    int n_edge;
    char tbuf[500];

    printf("XXXXX  \n");
    

        /*  printf("XXXXXXXXXXXXX \n"); */
/*     printf("%i %i %i \n", Total_N_markers, N_genomes, N_chromosomes); */
    fscanf(fp, "%s %i %*s %i", tbuf, &Total_N_markers, &N_genomes);
    if(N_genomes > 2){
        printf("Number of genomes is %i . Greater than 2 genomes not implemented.\n", N_genomes); //getchar();
    }
    { // read in first genome (0)
        i = 0;
        i_marker = 0;
        
        fscanf(fp, "%s %i", tempstr, &i_genome); // read in genome number
              printf("(Genome) tempstr, i_genome: [%s] [%i] \n", tempstr, i_genome);
        
        if(i != i_genome){printf("in get_data_from_file; i, i_genome: %i %i \n", i, i_genome);}
        
        fscanf(fp, "%s %i", temp2, &distances_this_genome);
            printf("(Distances) temp2: %s. distances_this_genome: %i \n", temp2, distances_this_genome);
        distances_in_file = distances_this_genome;
           
        fscanf(fp, "%s %i", tempstr, &N_chromosomes); // read in number of chromosomes for this genome
            //   printf("(Chromosomes)  %s %i \n", tempstr, N_chromosomes);
	printf("number of chromosomes in genome %i is %s %i \n", i_genome, tempstr, N_chromosomes);
      
        for(j=0; j<N_chromosomes; j++){
                //  printf("chromosome loop; j: %i \n", j);
            if(distances_this_genome == TRUE){
                    //    fscanf(fp, "%*s %i %lf", &i_chromosome, chromosome_length+j); // read number of this chromosome, and its length
                fscanf(fp, "%s %i %lf", tempstr, &i_chromosome, chromosome_length[i]+j); // read number of this chromosome, and its length 
                     printf("j, tempstr, chromosome#: %i  %s  %i. length: %g \n", j, tempstr, i_chromosome,chromosome_length[i][j]);
            }
            else{
                    //   chromosome_length[j] = N_markers[i][j]+1; 
	      fscanf(fp, "%s %i", tempstr, &i_chromosome); // read number of this chromosome
	      printf("AAAAAAAAAAAAAAAAA\n");
	      printf("tempstr chromosome#, j: [%s]   %i %i\n", tempstr, i_chromosome, j);
            } 
                // if(i_chromosome != j+1)
                // {printf("j, i_chromosome: %i %i \n", j, i_chromosome);}
            fscanf(fp, "%*s %i", &N_markers[i][j]); // read in number of markers for this chromosome
                   printf("chromosome#, j: %i %i, markers on this chrom: %i \n", i_chromosome, j, N_markers[i][j]);
            printf("number of markers on chromosome %i is %i \n", j, N_markers[i][j]);
            if(!distances_this_genome || !Use_distances) chromosome_length[i][j] = N_markers[i][j]+1; 
              
            fscanf(fp, "%s", tempstr); printf("::: [%s]\n", tempstr); // getchar();
            for(k=0; k<N_markers[i][j]; k++){
                char first[2] = {'\0', '\0'};
                char temp[32];
                if(TRUE){ // marker entries in  file are like:  -XK3   (i.e. no space between sign and rest of name. so works for signed numbers)
		  if(distances_this_genome) {
		    fscanf(fp, "%s %lf", temp, &dist);  // dist is a double
		    printf("ZZZZQQQ:   %s %g\n", temp, dist);
		    if(dist < 0.0 || dist > chromosome_length[i][j]){
		      printf("j, chromosome_length[j], dist: %i %g %g \n", j, chromosome_length[i][j], dist); getchar();}
		    if(!Use_distances) dist = k+1;
		    g[i][i_marker].marker_distance_on_chromosome = dist;
		    g[i][i_marker].marker_distance = cumulative_length_of_chromosomes + dist;
		    //  printf("j: %i chrom_length[j]: %g cume chrom length, dist: %g  %g \n",
		    //        j, chromosome_length[i][j], cumulative_length_of_chromosomes, dist); 
		  } else {
                        fscanf(fp, "%s", temp); dist = k+1;  // just use numbers 1,2,3... as distances from chrom end
			
                        g[i][i_marker].marker_distance_on_chromosome = dist;
                        g[i][i_marker].marker_distance = cumulative_length_of_chromosomes + dist;
                            //  printf("j: %i chrom_length[j]: %g cume chrom length, dist: %g  %g \n",
                            //      j, chromosome_length[i][j], cumulative_length_of_chromosomes, dist);
                    }
                   
                   
                    strncpy(first, temp, 1);
		    // fprintf(stderr, "ASDF: %s %s\n", first, temp);
                    if(strcmp(first, "-") == 0){ // sign is negative
                        strcpy(g[i][i_marker].marker_sign, "-");
                        strcpy(g[i][i_marker].marker_name, temp+1);
                            //  printf("%s %s\n", temp, g[i][i_marker].marker_name);
                    }
                    else if(strcmp(first, "+") == 0){ // sign is positive
                        strcpy(g[i][i_marker].marker_sign, "+");
                        strcpy(g[i][i_marker].marker_name, temp+1);
                    }
                    else if(strcmp(first, "?") == 0){ // sign is unknown
                        strcpy(g[i][i_marker].marker_sign, "?");
                        strcpy(g[i][i_marker].marker_name, temp+1);
                    }
                    else{ // no initial -,+, or ?. -> sign is positive, just copy all of temp to marker_name
                        strcpy(g[i][i_marker].marker_sign, "+");
                        strcpy(g[i][i_marker].marker_name, temp);
                    }
		    // fprintf(stderr, "SDFG: %i %i  %s %s \n", i, i_marker, g[i][i_marker].marker_sign, g[i][i_marker].marker_name); //getchar();
                }
                else{ // marker entries in  file are like:  - XK3
                    fscanf(fp, "%s", g[i][i_marker].marker_sign);
                    fscanf(fp, "%s", g[i][i_marker].marker_name);
                }
                    //  printf(" %s \n", g[i][i_marker].marker_name);
                g[i][i_marker].chromosome_number = j;
                g[i][i_marker].marker_position_on_chromosome = k;
                i_marker++;
            }//printf("\n");
            cumulative_length_of_chromosomes += chromosome_length[i][j];
                //   printf("*********j: %i cume length of chromosomes: %g . chromosome_length[j] %g\n",
                //      j, cumulative_length_of_chromosomes, chromosome_length[i][j]);
        } // end of loop over chromosomes

      /*   printf("first genome. L_g: %g \n", cumulative_length_of_chromosomes); */
/*         for(j=0; j<N_chromosomes; j++){ */
/*             printf("j, chrom_length[j]: %i  %g \n", j, chromosome_length[i][j]); */
/*         } */
        n_marker[i] = i_marker;
        qsort (g[i], i_marker, sizeof (Marker), marker_cmp);
        
        for(k=0; k<i_marker-1; k++){
            int compare;
            if((compare = strcmp(g[i][k].marker_name, g[i][k+1].marker_name)) >= 0){ // check that names are all distinct and in order
                printf("%s %s  %i \n", g[i][k].marker_name, g[i][k+1].marker_name, compare);
            }
        } printf("zxzxzxz \n"); //getchar();
        for(k=0; k<i_marker; k++){
            if(strcmp(g[i][k].marker_sign, "+") == 0) g[i][k].marker_number = k+1; // set marker_number according to lexicographic order of names
            else if(strcmp(g[i][k].marker_sign, "-") == 0) g[i][k].marker_number = -(k+1);
            else printf("marker sign in input file is neither + nor -. %s ./", g[i][k].marker_sign);
        }
           
    } // end read in genome 0
    L_g = cumulative_length_of_chromosomes;
    n_edge = N_chromosomes + n_marker[0];
    marker_spacing = (L_g/(double)n_edge);
    cumulative_length_of_chromosomes = 0.0;
    
      { // read in second genome (1)
          i = 1;
        i_marker = 0;
        
        fscanf(fp, "%s %i", tempstr, &i_genome); // read in genome number
        printf("(Genome) tempstr, i_genome: %s %i \n", tempstr, i_genome);
        
        if(i != i_genome){printf("in get_data_from_file; i, i_genome: %i %i \n", i, i_genome);}
        
        fscanf(fp, "%s %i", temp2, &distances_this_genome);
        printf("(Distances) temp2: %s. distances_this_genome: %i \n", temp2, distances_this_genome);
        if(distances_this_genome){ printf("Second genome in file has distances. Ignoring them. \n"); }
           
        fscanf(fp, "%s %i", tempstr, &N_chromosomes); // read in number of chromosomes for this genome
        printf("(Chromosomes)  %s %i \n", tempstr, N_chromosomes);
        printf("number of chromosomes in genome %i is %i \n", i_genome, N_chromosomes); 
      
        for(j=0; j<N_chromosomes; j++){
            printf("chromosome loop; j: %i \n", j);
            if(distances_this_genome == TRUE){
                    //    fscanf(fp, "%*s %i %lf", &i_chromosome, chromosome_length+j); // read number of this chromosome, and its length
                fscanf(fp, "%s %i %lf", tempstr, &i_chromosome, chromosome_length[i]+j); // read number of this chromosome, and its length 
                printf("j, tempstr, chromosome#: %i  %s  %i. length: %g \n", j, tempstr, i_chromosome,chromosome_length[i][j]);
            }
            else{
                    //   chromosome_length[j] = N_markers[i][j]+1; 
                fscanf(fp, "%*s %i", &i_chromosome); // read number of this chromosome
                    //  printf("chromosome#, j: %i %i, markers on this chrom: %i \n", i_chromosome, j, N_markers[i][j]);
            } 
                // if(i_chromosome != j+1)
           
            fscanf(fp, "%*s %i", &N_markers[i][j]); // read in number of markers for this chromosome
                //  if(i_chromosome != j+1)
            printf("number of markers on chromosome %i is %i \n", j, N_markers[i][j]);
            chromosome_length[i][j] = marker_spacing*(N_markers[i][j]+1);
            printf("marker_spacing: %g. chromosome_length[i][j]: %g \n", marker_spacing, chromosome_length[i][j]); 
            
            fscanf(fp, "%*s"); // read and discard 'Markers:'
            for(k=0; k<N_markers[i][j]; k++){
                char first[2] = {'\0', '\0'};
                char temp[32];
                if(TRUE){ // marker entries in  file are like:  -XK3   (i.e. no space between sign and rest of name. so works for signed numbers)
                    if(distances_this_genome) {
                        fscanf(fp, "%s %lf", temp, &dist);  // dist is a double
                        printf("%s %g\n", temp, dist);
                        if(dist < 0.0 || dist > chromosome_length[i][j]){
                            printf("j, chromosome_length[j], dist: %i %g %g \n", j, chromosome_length[i][j], dist); }
                    }
                    else {
                        fscanf(fp, "%s", temp); 
                    }
                        // read in dist, but don't use
                    dist = marker_spacing*(double)(k+1);  // just use numbers 1,2,3... as distances from chrom end
                    g[i][i_marker].marker_distance_on_chromosome = dist;
                    g[i][i_marker].marker_distance = cumulative_length_of_chromosomes + dist;
                        //   printf("j: %i chrom_length[j]: %g cume chrom length, dist: %g  %g \n",
                        //          j, chromosome_length[i][j], cumulative_length_of_chromosomes, dist);
                   
                   
                    strncpy(first, temp, 1);
                    if(strcmp(first, "-") == 0){ // sign is negative
                        strcpy(g[i][i_marker].marker_sign, "-");
                        strcpy(g[i][i_marker].marker_name, temp+1);
                            //  printf("%s %s\n", temp, g[i][i_marker].marker_name);
                    }
                    else if(strcmp(first, "+") == 0){ // sign is positive
                        strcpy(g[i][i_marker].marker_sign, "+");
                        strcpy(g[i][i_marker].marker_name, temp+1);
                    }
                    else if(strcmp(first, "?") == 0){ // sign is unknown
                        strcpy(g[i][i_marker].marker_sign, "?");
                        strcpy(g[i][i_marker].marker_name, temp+1);
                    }
                    else{ // no initial -,+, or ?. -> sign is positive, just copy all of temp to marker_name
                        strcpy(g[i][i_marker].marker_sign, "+");
                        strcpy(g[i][i_marker].marker_name, temp);
                    }
		    // fprintf(stderr, "SSSSDFG: %i %i  %s %s \n", i, i_marker, g[i][i_marker].marker_sign, g[i][i_marker].marker_name); //getchar();
                }
                else{ // marker entries in  file are like:  - XK3
                    fscanf(fp, "%s", g[i][i_marker].marker_sign);
                    fscanf(fp, "%s", g[i][i_marker].marker_name);
                }
                printf(" %s \n", g[i][i_marker].marker_name);
                g[i][i_marker].chromosome_number = j;
                g[i][i_marker].marker_position_on_chromosome = k;
                i_marker++;
            }printf("\n");
            cumulative_length_of_chromosomes += chromosome_length[i][j];
            printf("*********j: %i cume length of chromosomes: %g . chromosome_length[j] %g\n",
                   j, cumulative_length_of_chromosomes, chromosome_length[i][j]);
        } // end of loop over chromosomes

  printf("second genome. L_g: %g \n", cumulative_length_of_chromosomes);
        for(j=0; j<N_chromosomes; j++){
            printf("j, chrom_length[j]: %i  %g \n", j, chromosome_length[i][j]);
        }
        
        n_marker[i] = i_marker;
        qsort (g[i], i_marker, sizeof (Marker), marker_cmp);
        
        for(k=0; k<i_marker-1; k++){
            int compare;
            if((compare = strcmp(g[i][k].marker_name, g[i][k+1].marker_name)) >= 0){ // check that names are all distinct and in order
                printf("%s %s  %i \n", g[i][k].marker_name, g[i][k+1].marker_name, compare);
            }
        }
        for(k=0; k<i_marker; k++){
            if(strcmp(g[i][k].marker_sign, "+") == 0) g[i][k].marker_number = k+1; // set marker_number according to lexicographic order of names
            else if(strcmp(g[i][k].marker_sign, "-") == 0) g[i][k].marker_number = -(k+1);
            else printf("marker sign in input file is neither + nor -. %s ./", g[i][k].marker_sign);
        }
           
    } // end of read in genome 2
      
    if(n_marker[0] != n_marker[1]){
        printf("In get_data_from_file_new. Numbers of markers in two genomes are not the same: %i %i \n", n_marker[0], n_marker[1]);
        process_error();
    }
    if(!same_markers(n_marker[0], g[0], g[1])){
        printf("In get_data_from_file. Sets of markers in the two genomes are not the same. \n");
        process_error();
    }
    for(i=0; i<N_genomes; i++){
        for(k=0; k<i_marker; k++){ // get perm_array
            Marker* m;
            m = g[i]+k;     
                //   perm_array[i][m->chromosome_number][m->marker_position_on_chromosome] = m->marker_number;
            perm_array[i][m->chromosome_number][m->marker_position_on_chromosome] = *m;
        }
            //  print_perm_array(N_chromosomes, Nmarkers, perm_array[i]);
        the_perm = init_genome(N_chromosomes, N_markers[i], chromosome_length[i],
                               perm_array[i], cumulative_length_of_chromosomes);
        genomes[i] = the_perm;
        printf("genome %i : \n", i);
        print_perm(stdout, genomes[i]); 
    }
 /*    printf("bottom of get_data_from_file_new. genomes: \n"); */
/*      printf("first genome. L_g, A_i, A_t: %g %g %g \n", genomes[0]->L_g, genomes[0]->A_i, genomes[0]->A_t); */
/*      print_perm(stdout, genomes[0]); printf("end of first genome \n"); */
     
/*     printf("second genome. L_g, A_i, A_t: %g %g %g \n", genomes[1]->L_g, genomes[1]->A_i, genomes[1]->A_t); */
/*     print_perm(stdout, genomes[1]); printf("end of second genome \n"); */
    // getchar();
    return distances_in_file;
} // end of get_genomes_from_file_new




int
get_genome_from_file(FILE* fp, Permutation** genome, int Use_distances)
{
        // if Use_distances is TRUE and genome in file has distances
        // uses them, and spaces markers in second genome equally along genome, with
        // total length same as in first genome
        // else just use L_g = N+m for both, equally spaced (i.e. spacing is 1.0)
        // return value is TRUE if at least one genome has distance information
    int Total_N_markers, N_genomes, N_chromosomes;
    int N_markers[MAX_N_CHROMOSOMES] = {0}; // N_markers[i] = number of markers on chromosome i
    Marker perm_array[MAX_N_CHROMOSOMES][MAX_N_MARKERS_PER_CHROMOSOME];
    int j,k, i_genome, i_chromosome;
    Permutation* the_perm;
    Marker g[NMARKERMAX];
    int i_marker, n_marker = 0;
    int distances_in_file = FALSE, distances_this_genome; // TRUE/FALSE, indicates whether file has marker distances
    double dist;
    char tempstr[128], temp2[128];
    double chromosome_length[MAX_N_CHROMOSOMES];
    double cumulative_length_of_chromosomes = 0.0;
    double L_g, marker_spacing;
    int n_edge;
  
    fscanf(fp, "%*s %i %*s %i", &Total_N_markers, &N_genomes);
    if(N_genomes > 1){
        printf("Number of genomes is %i . Greater than 1 genomes not implemented.\n", N_genomes); //getchar();
    }
    { // read in the genome
        i_marker = 0;
        
        fscanf(fp, "%s %i", tempstr, &i_genome); // read in genome number
        if(0 != i_genome){printf("in get_data_from_file; i_genome: %i \n", i_genome);}
        fscanf(fp, "%s %i", temp2, &distances_this_genome);
        distances_in_file = distances_this_genome;    
        fscanf(fp, "%s %i", tempstr, &N_chromosomes); // read in number of chromosomes for this genome
            //  printf("number of chromosomes in genome %i is %i \n", i_genome, N_chromosomes);
      
        for(j=0; j<N_chromosomes; j++){
            if(distances_this_genome == TRUE){
                fscanf(fp, "%s %i %lf", tempstr, &i_chromosome, chromosome_length+j); // read number of this chromosome, and its length
            }
            else{
                fscanf(fp, "%*s %i", &i_chromosome); // read number of this chromosome
            } 
              
            fscanf(fp, "%*s %i", &N_markers[j]); // read in number of markers for this chromosome
                
                //  printf("number of markers on chromosome %i is %i \n", j, N_markers[j]);
            if(!distances_this_genome || !Use_distances) chromosome_length[j] = N_markers[j]+1; 
            
            fscanf(fp, "%*s"); // getchar();
            for(k=0; k<N_markers[j]; k++){
                char first[2] = {'\0', '\0'};
                char temp[32];
                if(TRUE){ // marker entries in  file are like:  -XK3   (i.e. no space between sign and rest of name. so works for signed numbers)
                    if(distances_this_genome) {
                        fscanf(fp, "%s %lf", temp, &dist);  // dist is a double
                            //  printf("%s %g\n", temp, dist);
                        if(dist < 0.0 || dist > chromosome_length[j]){
                            printf("j, chromosome_length[j], dist: %i %g %g \n", j, chromosome_length[j], dist); }
                        if(!Use_distances) dist = k+1;
                         g[i_marker].marker_distance_on_chromosome = dist;
                        g[i_marker].marker_distance = cumulative_length_of_chromosomes + dist;
                    }
                    else {
                        fscanf(fp, "%s", temp); dist = k+1;  // just use numbers 1,2,3... as distances from chrom end
                        g[i_marker].marker_distance_on_chromosome = dist;
                        g[i_marker].marker_distance = cumulative_length_of_chromosomes + dist;
                    }
                   
                   
                    strncpy(first, temp, 1);
                    if(strcmp(first, "-") == 0){ // sign is negative
                        strcpy(g[i_marker].marker_sign, "-");
                        strcpy(g[i_marker].marker_name, temp+1);
                    }
                    else if(strcmp(first, "+") == 0){ // sign is positive
                        strcpy(g[i_marker].marker_sign, "+");
                        strcpy(g[i_marker].marker_name, temp+1);
                    }
                    else if(strcmp(first, "?") == 0){ // sign is unknown
                        strcpy(g[i_marker].marker_sign, "?");
                        strcpy(g[i_marker].marker_name, temp+1);
                    }
                    else{ // no initial -,+, or ?. -> sign is positive, just copy all of temp to marker_name
                        strcpy(g[i_marker].marker_sign, "+");
                        strcpy(g[i_marker].marker_name, temp);
                    }
                }
                else{ // marker entries in  file are like:  - XK3
                    fscanf(fp, "%s", g[i_marker].marker_sign);
                    fscanf(fp, "%s", g[i_marker].marker_name);
                }
                g[i_marker].chromosome_number = j;
                g[i_marker].marker_position_on_chromosome = k;
                i_marker++;
            }
            cumulative_length_of_chromosomes += chromosome_length[j];
        } // end of loop over chromosomes

        n_marker = i_marker;
        
            // can sort by lexicographic order of marker names
            //  qsort(g, n_marker, sizeof (Marker), marker_cmp);    
            //     for(k=0; k<n_marker-1; k++){
            //        int compare;
            //        if((compare = strcmp(g[k].marker_name, g[k+1].marker_name)) >= 0){ // check that names are all distinct and in order
            //           printf("%s %s  %i \n", g[k].marker_name, g[k+1].marker_name, compare);
            //        }
            //     }
        
        for(k=0; k<n_marker; k++){
            if(strcmp(g[k].marker_sign, "+") == 0) g[k].marker_number = k+1; // set marker_number according to lexicographic order of names
            else if(strcmp(g[k].marker_sign, "-") == 0) g[k].marker_number = -(k+1);
            else {printf("marker sign in input file is neither + nor -. %s ./", g[k].marker_sign); getchar(); }
        }
    }// end read in genome 0
    
    L_g = cumulative_length_of_chromosomes;
    n_edge = N_chromosomes + n_marker;
    marker_spacing = (L_g/(double)n_edge);
    cumulative_length_of_chromosomes = 0.0; 
       
    for(k=0; k<n_marker; k++){ // get perm_array
        Marker* m;
        m = g+k;
        perm_array[m->chromosome_number][m->marker_position_on_chromosome] = *m;
    }
     
        //  print_perm_array(N_chromosomes, N_markers, perm_array);
      
    the_perm = init_genome(N_chromosomes, N_markers, chromosome_length,
                           perm_array, cumulative_length_of_chromosomes);
    *genome = the_perm;
  
    return distances_in_file;
} // end of get_genome_from_file 


void
print_perm_array(int Nchrom, int* Nmarkers, Marker perm_array[][MAX_N_MARKERS_PER_CHROMOSOME])
{
    Marker* m;
    int i, j;
    for(i=0; i<Nchrom; i++){
        for(j=0; j<Nmarkers[i]; j++){
            m = perm_array[i]+j;
                //  printf("%4i %4i %g \n", m->marker_number, m->marker_position_on_chromosome, m->marker_distance_on_chromosome);
            printf("%5i %5i %5i  %s\n", m->marker_number, m->chromosome_number, m->marker_position_on_chromosome, m->marker_name);
        }printf("\n");
    }printf("\n");
}  // end of print_perm_array

/* int */
/* marker_cmp(const Marker* m1, const Marker* m2) */
/* { */
/*     return strcmp(m1->marker_name, m2->marker_name); */
/* } */

int
marker_cmp(const void* am1, const void* am2)
{
    Marker* m1 = (Marker*)am1, * m2 = (Marker*)am2;
    return strcmp(m1->marker_name, m2->marker_name);
}


int
same_markers(int N, Marker* g1, Marker* g2)
{
    int i;
   
    for(i=0; i<N; i++){
        printf("%s %s \n", g1[i].marker_name, g2[i].marker_name);
        if(strcmp(g1[i].marker_name, g2[i].marker_name) != 0){
            printf("In same_markers. Marker names dont match: %s %s \n", g1[i].marker_name, g2[i].marker_name);
             return FALSE;
        }
    }
    return TRUE;
} // end same_markers



double
get_delta_d_prob_new(Cycle_decomposition* the_cd, int delta_d, int same_eds, int is_inversion,
                     int n_inv, int n_trans, double Epsilon, double lambda_ratio, int dnsw)
{
        // Given a cycle decomposition, returns the probability of choosing
        // an particular inversion/translocation with the specified delta_d
   
    double prob_ddm1, prob_dd0, prob_dd1; //, prob_tot;
    double retval = -1.0;
    double p_inv, p_trans, denom;

    prob_ddm1 = step_fcn(the_cd->n_dd[0])*(1.0 - Epsilon);
        //  prob_dd0, prob_dd1;
        // prob_tot = prob_ddm1 + prob_dd0 + prob_dd1;
    
    if(delta_d == -1){
            //   printf("in get_delta_c_prob. prob_ddm1:
            //         printf("in get_delta_d_prob. delta_d: %i . n_ddm1: %i . prob_ddm1: %g \n", delta_d, the_cd->n_dd[0], prob_ddm1);
        if(ddm1_prefer_inversions == FALSE){
            if(use_dnswitches == FALSE){
                
                retval = prob_ddm1/(double)the_cd->n_dd[0];
                    //  printf("in get_delta_d_prob_new. retval, prob_ddm1, n_ddm1: %g %g %i\n", retval, prob_ddm1, the_cd->n_dd[0]);
            }
            else{
                double p[5];
                int i; //, dsw_flipped;
                    //  int dsw;
                int dnsw_index = dnsw+2;
                    //  int dsw1;
                 
                denom = 0.0;
                  
                for(i=0; i<5; i++){
                        //  printf("%i  ", (the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]));
                    if(STEPFCN_USEDNSW){
                        p[i] = dswfactor[i]*step_fcn(the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                    }
                    else{
                            //  p[i] = dswfactor[i]*(double)(the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                        p[i] = p_dnsw(i, the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                    }
                    denom += p[i];
                } for(i=0; i<5; i++){ p[i] /= denom; }
                    //  printf("\n");
              
                retval = prob_ddm1*p[dnsw_index]/(double)(the_cd->n_dsw_dcgm1[dnsw_index] + the_cd->n_dsw_dci1[dnsw_index]);
                    //     printf("in get_delta_d_prob_new: retval, prob_ddm1, dsw_index, p[dsw_index], denom: %g %g %i %g  %i \n",
                    //          retval, prob_ddm1, dnsw_index, p[dnsw_index], (the_cd->n_dsw_dcgm1[dnsw_index] + the_cd->n_dsw_dci1[dnsw_index]));
                  if(p[dnsw_index] == 0.0){
                      Cycle_decomposition acd;

                      
                      get_CD(the_cd->perm1, the_cd->perm2, &acd);
                      for(i=0; i<5; i++){
                          printf("%i  ", (acd.n_dsw_inv_ddm1[i] + acd.n_dsw_trans_ddm1[i])); 
                      }
                      printf(" \n");
                     
                      print_perm(stdout, the_cd->perm1); print_perm(stdout, the_cd->perm2);
                      getchar();
                  }
            }
        }
        else{ // prefer inversions to translations
            p_inv = DDM1_IT_PROB_RATIO*step_fcn(the_cd->n_dd_inv[0]);
            p_trans = step_fcn(the_cd->n_dd_trans[0]);
            denom = p_inv + p_trans;
            p_inv /= denom; p_trans /= denom;
            if(use_dnswitches == FALSE){
                if(is_inversion) retval = prob_ddm1*p_inv/(double)the_cd->n_dd_inv[0];
                else retval = prob_ddm1*p_trans/(double)the_cd->n_dd_trans[0];
            }
            else{
                double p[5];
                int i; //, dsw_flipped;
                    //  int dsw, dsw1;
                int dnsw_index;
                denom = 0.0;

                   
                dnsw_index = dnsw + 2;
                if(is_inversion) { // inversion
                    for(i=0; i<5; i++){
                        if(STEPFCN_USEDNSW){
                            p[i] = dswfactor[i]*step_fcn(the_cd->n_dsw_inv_ddm1[i]);
                        }
                        else{
                                //  p[i] = dswfactor[i]*(double)(the_cd->n_dsw_inv_ddm1[i]);
                            p[i] = p_dnsw(i, the_cd->n_dsw_inv_ddm1[i]);
                        }
                        denom += p[i];
                    } for(i=0; i<5; i++){ p[i] /= denom; }
                       /*      printf("inv: prob_ddm1, p_inv, dsw_index, p[dsw_index], the_cd->n_dsw_inv_ddm1[dsw_index]: %g %g %i %g %i \n", */
/*                            prob_ddm1, p_inv, dnsw_index, p[dnsw_index], the_cd->n_dsw_inv_ddm1[dnsw_index]); */
                    retval = prob_ddm1*p_inv*p[dnsw_index]/(double)the_cd->n_dsw_inv_ddm1[dnsw_index];
                }
                else { // translocation
                    for(i=0; i<5; i++){
                        if(STEPFCN_USEDNSW){
                            p[i] = dswfactor[i]*step_fcn(the_cd->n_dsw_trans_ddm1[i]);
                        }
                        else{
                                //  p[i] = dswfactor[i]*(double)(the_cd->n_dsw_trans_ddm1[i]);
                           p[i] = p_dnsw(i, the_cd->n_dsw_trans_ddm1[i]); 
                        }
                        denom += p[i];
                    } for(i=0; i<5; i++){ p[i] /= denom; }
                  /*    printf("trans: prob_ddm1, p_trans, dsw_index, p[dsw_index], the_cd->n_dsw_inv_ddm1[dsw_index]: %g %g %i %g %i \n", */
/*                             prob_ddm1, p_trans, dnsw_index, p[dnsw_index], the_cd->n_dsw_inv_ddm1[dnsw_index]);  */
                    retval = prob_ddm1*p_trans*p[dnsw_index]/(double)the_cd->n_dsw_trans_ddm1[dnsw_index];
                }
            }
        }   
    }
    else{ // delta_d = 0 or 1

        the_cd->n_dcg[2] = get_n_dcg1(the_cd); // also does the_cd->n00_deds1       
        the_cd->n00_deds2 = get_n00_disteds2(the_cd);
            
        the_cd->n_dd[1] = the_cd->n00_seds + the_cd->n00_deds1 + the_cd->n00_deds2;
        the_cd->n_dd_inv[1] = the_cd->n00_seds_inv + the_cd->n00_deds1_inv + the_cd->n00_deds2_inv;
        the_cd->n_dd_trans[1] = the_cd->n00_seds_trans + the_cd->n00_deds1_trans + the_cd->n00_deds2_trans;
            
        the_cd->n_dci[0] = n_inv + n_trans - (the_cd->n_dcg[0] + the_cd->n_dcg[2] + the_cd->n_dci[2] + the_cd->n_dd[1]);
        the_cd->n_dd[2] = the_cd->n_dcg[2] + the_cd->n_dci[0];
            
        the_cd->n_dci_inv[0] = n_inv - (the_cd->n_dcg_inv[0] + the_cd->n_dcg_inv[2] + the_cd->n_dci_inv[2] + the_cd->n_dd_inv[1]);
        the_cd->n_dd_inv[2] = the_cd->n_dcg_inv[2] + the_cd->n_dci_inv[0];
            
        the_cd->n_dci_trans[0] = n_trans - (the_cd->n_dcg_trans[0] + the_cd->n_dcg_trans[2] + the_cd->n_dci_trans[2] + the_cd->n_dd_trans[1]);
        the_cd->n_dd_trans[2] = the_cd->n_dcg_trans[2] + the_cd->n_dci_trans[0];

        prob_dd0 = 1.0 - prob_ddm1;
        prob_dd1 = step_fcn(the_cd->n_dd[2])*prob_dd0*EPSILON2FACTOR;
        prob_dd0 /= 1.0 + step_fcn(the_cd->n_dd[2])*EPSILON2FACTOR;
        prob_dd1 /= 1.0 + step_fcn(the_cd->n_dd[2])*EPSILON2FACTOR;
    
        if(delta_d == 1){
               /*       printf("in get_delta_d_prob. delta_d: %i . n_dds: %i %i %i. prob_dd0/1: %g %g. n00s: %i %i %i \n", */
/*                             delta_d, the_cd->n_dd[0], the_cd->n_dd[1], the_cd->n_dd[2], prob_dd0, prob_dd1, */
/*                             the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
            if(PREFER_INVERSIONS == TRUE){ // prefer inversions
                if(USE_IPF_LAMBDA_RATIO){
                    p_inv = pow(lambda_ratio, IPF_EXPONENT)*the_cd->n_dd_inv[2];
                }
                else{p_inv = pow(1.0, IPF_EXPONENT)*the_cd->n_dd_inv[2]; }
                p_trans = 1.0*the_cd->n_dd_trans[2];
                denom = p_inv + p_trans;
                p_inv /= denom; p_trans /= denom; // normalize
                if(is_inversion){
                    retval = prob_dd1*p_inv/(double)the_cd->n_dd_inv[2];
                }
                else{
                    retval = prob_dd1*p_trans/(double)the_cd->n_dd_trans[2];
                }
                    //  printf("in delta_d_prob_new. delta d: %i  is_inversion: %i  p_inv,trans: %g %g  retval: %g \n", delta_d, is_inversion, p_inv, p_trans, retval);
            }
            else{
                retval =  prob_dd1/(double)the_cd->n_dd[2];
            }
        }
        else if(delta_d == 0){
            if(PREFERDD0SEDS && (the_cd->n_dd[0] <= 0)){
                double p_dd0seds = (the_cd->n00_seds > 0)? PDD0SEDS: 0;
                double p_dd0deds = 1.0 - p_dd0seds;
            
                if(same_eds){
                        /*   printf("  In get_delta_d_prob. delta_d=0, seds, prob_dd0, prob_tot, n00_seds,deds12 %g %g %i %i %i \n", */
/*                        prob_dd0, prob_tot, the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
                    retval = (prob_dd0)*p_dd0seds/(double)the_cd->n00_seds;
                }
                else{
                        /*   printf("  In get_delta_d_prob. delta_d=0, deds, prob_dd0/prob_tot, p_dd-deds, n00_seds,deds12 %g %g %i %i %i\n", */
/*                        prob_dd0/prob_tot, p_dd0deds, the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
                    retval = (prob_dd0)*p_dd0deds/(double)(the_cd->n00_deds1 + the_cd->n00_deds2); 
                }
            }
            else{
                if(PREFER_INVERSIONS == TRUE){ // prefer inversions
                    if(USE_IPF_LAMBDA_RATIO){
                        p_inv = pow(lambda_ratio, IPF_EXPONENT)*the_cd->n_dd_inv[1];
                    } 
                    else{ p_inv = pow(1.0, IPF_EXPONENT)*the_cd->n_dd_inv[1]; }
                    p_trans = 1.0*the_cd->n_dd_trans[1];
                    denom = p_inv + p_trans;
                    p_inv /= denom; p_trans /= denom; // normalize
                    if(is_inversion){
                        retval = prob_dd0*p_inv/(double)the_cd->n_dd_inv[1];
                    }
                    else{
                        retval = prob_dd0*p_trans/(double)the_cd->n_dd_trans[1];
                    }
                    
                }
                else{
                    retval = prob_dd0/(double)the_cd->n_dd[1];
                }
            }
            /*  printf("in get_delta_d_prob. delta_d: %i . n_dds: %i %i %i. prob_dd0/1: %g %g. n00s: %i %i %i \n", */
/*                             delta_d, the_cd->n_dd[0], the_cd->n_dd[1], the_cd->n_dd[2], prob_dd0, prob_dd1, */
/*                             the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
        }
        else
            printf("in get_delta_d_prob. delta_d: %i . Should be -1, 0 or 1.\n", delta_d);
    }
    
        //  printf("in get_delta_d_prob. delta_d, prob: %i %g \n", delta_d, retval);
    if(retval < 0.0) {printf("in get_delta_d_prob_new. retval < 0: %g \n", retval); }
    return retval;
} // end of get_delta_d_prob_new

Reversal
delta_ci_cg_reverse_new(Cycle_decomposition* the_cd, double epsilon, int n_inv, int n_trans, double* prob, double lambda_ratio)
{
#define PRINTSTUFF FALSE
        // given a cycle decomposition, and epsilon.
        // return reversal
    
    int n_dcgm1 = the_cd->n_dcg[0],
        n_dci1 = the_cd->n_dci[2],
        n_dd0, n_dcg1, n_dcim1;
    int n;
    double p_ddm1 = step_fcn(n_dcgm1 + n_dci1)*(1.0 - epsilon);
    double p_ddge0 = 1.0 - p_ddm1;
    double p_dd0, p_dd1;
    double randnumb = drand(&rng_long);
    Reversal rev;
    int distance = the_cd->n_mark - the_cd->n_chrom - the_cd->n_int_cycles + the_cd->c_g;
    double p_00seds, p_00deds;
    double p_inv, p_trans;
    double denom;

  /*   printf("%i %i \n", n_dcgm1, n_dci1); */
/*     printf("epsilon, p_ddm1, p_ddge0:  %g  %g  %g  \n", epsilon, p_ddm1, p_ddge0); */
        //   printf("d: %i  n_ddm1, inv, trans, total: %i %i %i \n", the_cd->d, the_cd->n_dd_inv[0], the_cd->n_dd_trans[0], the_cd->n_dd[0]);
    if(PRINTSTUFF) printf("%12g  %4i ", epsilon, the_cd->n_mark - the_cd->n_chrom - the_cd->n_int_cycles + the_cd->c_g);
    if(PRINTSTUFF) printf("%5i %5i %5i %5i %5i ", n_dci1, n_dcgm1, n_dd0, n_dcim1, n_dcg1);
    if(PRINTSTUFF) printf("%12g %12g %12g  ", p_ddm1, p_dd0, p_dd1);
    if(randnumb < p_ddm1){ // usual "downhill step" branch, i.e. delta_d = -1
        if(ddm1_prefer_inversions == FALSE){
            if(use_dnswitches == FALSE){
                n = (int)(drand(&rng_long)*(n_dcgm1 + n_dci1));
                if(n < n_dcgm1){ // do a delta c_g = -1 step
                    rev = get_spec_dcgm1_rev(the_cd, n+1, EITHER);
                    rev.delta_c_i = 0; rev.delta_c_g = -1; rev.same_eds = FALSE;
                }
                else{ // do a delta c_i = +1 step
                    n -= n_dcgm1;
                    rev = get_spec_dci1_rev(the_cd, n+1, EITHER);
                        //   rev = get_spec_dci1_rev_old(the_cd, n+1);
                    rev.delta_c_i = 1; rev.delta_c_g = 0; rev.same_eds = TRUE; 
                }
                *prob = p_ddm1/(double)(n_dcgm1 + n_dci1);
                    //  printf("in get_ci_cg_revse_new. *prob, prob_ddm1, n_ddm1: %g %g %i\n", *prob, p_ddm1, n_dcgm1+n_dci1);
            }
            else { // prefer inv/trans which reduce n_switches
                int i;
                double p[5], the_p;
                randnumb = drand(&rng_long);
                denom = 0.0;
                for(i=0; i<5; i++){
                    if(STEPFCN_USEDNSW){
                    p[i] = dswfactor[i]*step_fcn(the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                    }
                    else {
                            //   p[i] = dswfactor[i]*(double)(the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                        p[i] = p_dnsw(i, the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                    }
                        //   printf(" %i " , the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                    denom += p[i];
                } // printf("\n");
             /*    for(i=0; i<5; i++){ printf(" %i ", the_cd->n_dsw_dcgm1[i]); } printf("\n"); */
/*                  for(i=0; i<5; i++){ printf(" %i ", the_cd->n_dsw_dci1[i]); } printf("\n"); */
                for(i=0; i<5; i++){ p[i] /= denom; }
                for(i=0, the_p =0.0; i<5; i++){
                    the_p += p[i];
                    if(randnumb < the_p){ // do a delta_n_sw = i-2 inv/trans
                        n = (int)(drand(&rng_long)*(double)(the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i]));
                        if(n < the_cd->n_dsw_dcgm1[i]){ // dcgm1 step
                                //     printf("in get_delta_ci... before get_spec_dcgm1_rev_new. dsw_index: %i   n: %i \n", i, n);
                            rev = get_spec_dcgm1_rev_new(the_cd, n+1, EITHER, i);
                                //    printf("in get_delta_ci... after get_spec_dcgm1_rev_new. \n"); 
                        }
                        else{
                                //   printf(" n, n_dsw_dcgm1[i], i: %i %i %i   ", n, the_cd->n_dsw_dcgm1[i], i);
                            n -= the_cd->n_dsw_dcgm1[i];
                                //   printf(" n: %i  n_dsw_dci1[i]: %i \n", n, the_cd->n_dsw_dci1[i]);
                            assert(n < the_cd->n_dsw_dci1[i]);
                                //    printf("in get_delta_ci... before get_spec_dci1_rev_new1. dsw_index: %i   n: %i \n", i, n);
                            rev = get_spec_dci1_rev_new1(the_cd, n+1, EITHER, i);
                                //    printf("in get_delta_ci... after get_spec_dci1_rev_new1. \n");
                        }
                    /*     printf("in get_delta_ci_cg... i, prob_ddm1, p[i], (the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i]): %i %g %g  %i \n", */
/*                                i, p_ddm1, p[i], (the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i])); */
                        *prob = p_ddm1*p[i]/(double)(the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i]);
                        i = 5; // so this is last time through loop
                    }
                }
              
                
            }
        }
        else { // prefer inversions over translocations
            p_inv = DDM1_IT_PROB_RATIO*step_fcn(the_cd->n_dd_inv[0]);
            p_trans = step_fcn(the_cd->n_dd_trans[0]);
            denom = p_inv + p_trans;
            p_inv /= denom; p_trans /= denom;
            if(use_dnswitches == FALSE){
                if(drand(&rng_long) < p_inv){ // do an inversion with delta_d=-1
                    n = (int)(drand(&rng_long)*(the_cd->n_dd_inv[0]));
                    if(n < the_cd->n_dcg_inv[0]){ // do a delta c_g = -1 step            
                                                
                        rev = get_spec_dcgm1_rev(the_cd, n+1, INVERSION);
                        rev.delta_c_i = 0; rev.delta_c_g = -1; rev.same_eds = FALSE;
                    }
                    else{ // do a delta c_i = +1 step
                        n -= the_cd->n_dcg_inv[0];
                            //  printf("n, the_cd->n_dci_inv[2]: %i %i \n", n, the_cd->n_dci_inv[2]);
                        rev = get_spec_dci1_rev(the_cd, n+1, INVERSION);
                        rev.delta_c_i = 1; rev.delta_c_g = 0; rev.same_eds = TRUE; 
                    }
                    *prob = p_ddm1*p_inv/(double)the_cd->n_dd_inv[0];
                }
                else{ // do a translocation with delta_d=-1
                    n = (int)(drand(&rng_long)*(the_cd->n_dd_trans[0]));
                    if(n < the_cd->n_dcg_trans[0]){ // do a delta c_g = -1 step
                        rev = get_spec_dcgm1_rev(the_cd, n+1, TRANSLOCATION);
                        rev.delta_c_i = 0; rev.delta_c_g = -1; rev.same_eds = FALSE;
                    }
                    else{ // do a delta c_i = +1 step
                        n -= the_cd->n_dcg_trans[0];
                            //  printf("n, the_cd->n_dci_trans[2]: %i %i \n", n, the_cd->n_dci_trans[2]);
                        rev = get_spec_dci1_rev(the_cd, n+1, TRANSLOCATION);
                        rev.delta_c_i = 1; rev.delta_c_g = 0; rev.same_eds = TRUE; 
                    }
                    *prob = p_ddm1*p_trans/(double)the_cd->n_dd_trans[0];
                }
            }
            else{ // use dnsw
              
                if(drand(&rng_long) < p_inv){ // do an inversion with delta_d=-1

                    int i;
                    double p[5], the_p;
                    randnumb = drand(&rng_long);
                    denom = 0.0;
                    for(i=0; i<5; i++){
                        if(STEPFCN_USEDNSW){
                            p[i] = dswfactor[i]*step_fcn(the_cd->n_dsw_inv_ddm1[i]);
                        }
                        else {
                                //  p[i] = dswfactor[i]*(double)(the_cd->n_dsw_inv_ddm1[i]);
                            p[i] = p_dnsw(i, the_cd->n_dsw_inv_ddm1[i]); 
                        }
                            //   printf(" %i " , the_cd->n_dsw_inv_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                        denom += p[i];
                    } // printf("\n");
                        /*    for(i=0; i<5; i++){ printf(" %i ", the_cd->n_dsw_dcgm1[i]); } printf("\n"); */
/*                  for(i=0; i<5; i++){ printf(" %i ", the_cd->n_dsw_dci1[i]); } printf("\n"); */
                    for(i=0; i<5; i++){ p[i] /= denom; }
                    for(i=0, the_p =0.0; i<5; i++){
                        the_p += p[i];
                        if(randnumb < the_p){ // do a delta_n_sw = i-2 inv
                            n = (int)(drand(&rng_long)*(double)(the_cd->n_dsw_inv_dcgm1[i] + the_cd->n_dsw_inv_dci1[i]));
                            if(n < the_cd->n_dsw_inv_dcgm1[i]){ // dcgm1 step
                                    //     printf("in get_delta_ci... before get_spec_dcgm1_rev_new. dsw_index: %i   n: %i \n", i, n);
                                rev = get_spec_dcgm1_rev_new(the_cd, n+1, INVERSION, i);
                                    //    printf("in get_delta_ci... after get_spec_dcgm1_rev_new. \n"); 
                            }
                            else{
                                    //   printf(" n, n_dsw_dcgm1[i], i: %i %i %i   ", n, the_cd->n_dsw_dcgm1[i], i);
                                n -= the_cd->n_dsw_inv_dcgm1[i];
                                    //   printf(" n: %i  n_dsw_dci1[i]: %i \n", n, the_cd->n_dsw_dci1[i]);
                                assert(n < the_cd->n_dsw_inv_dci1[i]);
                                    //    printf("in get_delta_ci... before get_spec_dci1_rev_new1. dsw_index: %i   n: %i \n", i, n);
                                rev = get_spec_dci1_rev_new1(the_cd, n+1, INVERSION, i);
                                    //    printf("in get_delta_ci... after get_spec_dci1_rev_new1. \n");
                            }
                                /*     printf("in get_delta_ci_cg... i, prob_ddm1, p[i], (the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i]): %i %g %g  %i \n", */
/*                                i, p_ddm1, p[i], (the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i])); */
                            *prob = p_ddm1*p_inv*p[i]/(double)(the_cd->n_dsw_inv_dcgm1[i] + the_cd->n_dsw_inv_dci1[i]);
                            i = 5; // so this is last time through loop
                        }
                    }
                      
                }
                else{ // do a translocation with delta_d=-1
                   
                    int i;
                    double p[5], the_p;
                    randnumb = drand(&rng_long);
                    denom = 0.0;
                    for(i=0; i<5; i++){
                        if(STEPFCN_USEDNSW){
                            p[i] = dswfactor[i]*step_fcn(the_cd->n_dsw_trans_ddm1[i]);
                        }
                        else {
                                //  p[i] = dswfactor[i]*(double)(the_cd->n_dsw_trans_ddm1[i]);
                            p[i] = p_dnsw(i, the_cd->n_dsw_trans_ddm1[i]);
                        }
                            //   printf(" %i " , the_cd->n_dsw_trans_ddm1[i] + the_cd->n_dsw_trans_ddm1[i]);
                        denom += p[i];
                    } // printf("\n");
                        /*    for(i=0; i<5; i++){ printf(" %i ", the_cd->n_dsw_dcgm1[i]); } printf("\n"); */
/*                  for(i=0; i<5; i++){ printf(" %i ", the_cd->n_dsw_dci1[i]); } printf("\n"); */
                    for(i=0; i<5; i++){ p[i] /= denom; }
                    for(i=0, the_p =0.0; i<5; i++){
                        the_p += p[i];
                        if(randnumb < the_p){ // do a delta_n_sw = i-2 trans
                            n = (int)(drand(&rng_long)*(double)(the_cd->n_dsw_trans_dcgm1[i] + the_cd->n_dsw_trans_dci1[i]));
                            if(n < the_cd->n_dsw_trans_dcgm1[i]){ // dcgm1 step
                                    //     printf("in get_delta_ci... before get_spec_dcgm1_rev_new. dsw_index: %i   n: %i \n", i, n);
                                rev = get_spec_dcgm1_rev_new(the_cd, n+1, TRANSLOCATION, i);
                                    //    printf("in get_delta_ci... after get_spec_dcgm1_rev_new. \n"); 
                            }
                            else{
                                    //   printf(" n, n_dsw_dcgm1[i], i: %i %i %i   ", n, the_cd->n_dsw_dcgm1[i], i);
                                n -= the_cd->n_dsw_trans_dcgm1[i];
                                    //   printf(" n: %i  n_dsw_dci1[i]: %i \n", n, the_cd->n_dsw_dci1[i]);
                                assert(n < the_cd->n_dsw_trans_dci1[i]);
                                    //    printf("in get_delta_ci... before get_spec_dci1_rev_new1. dsw_index: %i   n: %i \n", i, n);
                                rev = get_spec_dci1_rev_new1(the_cd, n+1, TRANSLOCATION, i);
                                    //    printf("in get_delta_ci... after get_spec_dci1_rev_new1. \n");
                            }
                                /*     printf("in get_delta_ci_cg... i, prob_ddm1, p[i], (the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i]): %i %g %g  %i \n", */
/*                                i, p_ddm1, p[i], (the_cd->n_dsw_dcgm1[i] + the_cd->n_dsw_dci1[i])); */
                            *prob = p_ddm1*p_trans*p[i]/(double)(the_cd->n_dsw_trans_dcgm1[i] + the_cd->n_dsw_trans_dci1[i]);
                            i = 5; // so this is last time through loop
                        }
                    } 
                }
            }       
        
        } // end of prefer inversions over translocations
    } // end of delta d = -1 step branch
    else { // take a delta_d = 0 or +1 step
          
        the_cd->n_dcg[2] = get_n_dcg1(the_cd); // also does the_cd->n00_deds1       
        the_cd->n00_deds2 = get_n00_disteds2(the_cd);
            
        the_cd->n_dd[1] = the_cd->n00_seds + the_cd->n00_deds1 + the_cd->n00_deds2;
        the_cd->n_dd_inv[1] = the_cd->n00_seds_inv + the_cd->n00_deds1_inv + the_cd->n00_deds2_inv;
        the_cd->n_dd_trans[1] = the_cd->n00_seds_trans + the_cd->n00_deds1_trans + the_cd->n00_deds2_trans;
            
        the_cd->n_dci[0] = n_inv + n_trans - (the_cd->n_dcg[0] + the_cd->n_dcg[2] + the_cd->n_dci[2] + the_cd->n_dd[1]);
        the_cd->n_dd[2] = the_cd->n_dcg[2] + the_cd->n_dci[0];
            
        the_cd->n_dci_inv[0] = n_inv - (the_cd->n_dcg_inv[0] + the_cd->n_dcg_inv[2] + the_cd->n_dci_inv[2] + the_cd->n_dd_inv[1]);
        the_cd->n_dd_inv[2] = the_cd->n_dcg_inv[2] + the_cd->n_dci_inv[0];
            
        the_cd->n_dci_trans[0] = n_trans - (the_cd->n_dcg_trans[0] + the_cd->n_dcg_trans[2] + the_cd->n_dci_trans[2] + the_cd->n_dd_trans[1]);
        the_cd->n_dd_trans[2] = the_cd->n_dcg_trans[2] + the_cd->n_dci_trans[0];
                

        p_dd0 = p_ddge0/(1.0+step_fcn(the_cd->n_dd[2])*EPSILON2FACTOR);
        p_dd1 = p_ddge0*step_fcn(the_cd->n_dd[2])*EPSILON2FACTOR/(1.0+step_fcn(the_cd->n_dd[2])*EPSILON2FACTOR);
         
        n_dd0 = the_cd->n_dd[1];
        n_dcg1 = the_cd->n_dcg[2];
        n_dcim1 = the_cd->n_dci[0];
       
        
        if(randnumb < p_ddm1 + p_dd0){ // take a delta d = 0 step
            if(PREFERDD0SEDS && (the_cd->n_dd[0] <= 0)){ // getting close; hurdles, etc. may be problem, so prefer seds moves
                    // PDD0SEDS controls preference for seds steps
                p_00seds = (the_cd->n00_seds > 0)? PDD0SEDS: 0.0;
                p_00deds = 1.0 - p_00seds;
                if(drand(&rng_long) < p_00seds){ // do a 00_seds step
                    n = (int)(drand(&rng_long)*the_cd->n00_seds);
                        /*   printf("In delta_ci_cg. 00seds, p_dd0, p_00seds,deds12, nseds: %g %g %i %i %i \n", */
/*                        p_dd0, p_00seds, the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
                    rev = get_spec_00seds_rev(the_cd, n+1, EITHER);
                    *prob = p_dd0*p_00seds/(double)the_cd->n00_seds;
                    rev.same_eds = TRUE;
                }
                else{ // do a 00_deds step
                    n = (int)(drand(&rng_long)*(the_cd->n00_deds1 + the_cd->n00_deds2));
                        /*     printf("In delta_ci_cg. 00deds, p_dd0, p_00deds, nseds,deds1,2: %g %g %i %i %i \n", */
/*                        p_dd0, p_00deds, the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
                    if(n < the_cd->n00_deds1){
                        rev = get_spec_00deds1_rev(the_cd, n+1, EITHER);
                    }
                    else if((n -= the_cd->n00_deds1) < the_cd->n00_deds2){
                        rev = get_spec_00deds2_rev(the_cd, n+1, EITHER);
                    }
                    else printf("in delta_ci_cg_reverse. d = 0. Problem.\n");
                    *prob = p_dd0*p_00deds/(double)(the_cd->n00_deds1 + the_cd->n00_deds2);
                    rev.same_eds = FALSE;
                   /*  printf("delta_d=0, deds, p_dd0: %g p_00deds: %g  *prob: %g. n_00's: %i %i %i \n", */
/*                            p_dd0, p_00deds, *prob, the_cd->n00_seds, the_cd->n00_deds1, the_cd->n00_deds2); */
                }
            }
            else{ // don't prefer seds
                if(PREFER_INVERSIONS == TRUE){ // prefer inversions
                        //    double denom;
                    if(USE_IPF_LAMBDA_RATIO){
                        p_inv = pow(lambda_ratio, IPF_EXPONENT)*the_cd->n_dd_inv[1];
                    }
                    else{ p_inv = pow(1.0, IPF_EXPONENT)*the_cd->n_dd_inv[1]; }
                        //     printf("IPF^IPF_exponent: %9g ", pow(INVERSION_PREFER_FACTOR, IPF_EXPONENT));
                    p_trans = 1.0*the_cd->n_dd_trans[1];
                    denom = p_inv + p_trans;
                    p_inv /= denom; p_trans /= denom; // normalize
                   
                    if(drand(&rng_long) < p_inv){ // choose from among the inversions
                        n = (int)(drand(&rng_long)*the_cd->n_dd_inv[1]);
                        assert(the_cd->n_dd_inv[1] == (the_cd->n00_seds_inv + the_cd->n00_deds1_inv + the_cd->n00_deds2_inv));
                        if(n < the_cd->n00_seds_inv){
                            rev = get_spec_00seds_rev(the_cd, n+1, INVERSION);
                            rev.same_eds = TRUE;
                        }
                        else if ((n -= the_cd->n00_seds_inv) < the_cd->n00_deds1_inv){
                            rev = get_spec_00deds1_rev(the_cd, n+1, INVERSION);
                            rev.same_eds = FALSE;
                        }
                        else if((n -= the_cd->n00_deds1_inv) < the_cd->n00_deds2_inv){
                            rev = get_spec_00deds2_rev(the_cd, n+1, INVERSION);
                            rev.same_eds = FALSE;
                        }
                        else printf("in delta_ci_cg_reverse. d = 0. Problem.\n");   
                            //    printf("dd=0, inversion. \n"); 
                        *prob = p_dd0*p_inv/(double)the_cd->n_dd_inv[1];
                    }
                    else{ // choose from among the translocations
                        n = (int)(drand(&rng_long)*the_cd->n_dd_trans[1]);
                        assert(the_cd->n_dd_trans[1] == (the_cd->n00_seds_trans + the_cd->n00_deds1_trans + the_cd->n00_deds2_trans));
                        if(n < the_cd->n00_seds_trans){
                            rev = get_spec_00seds_rev(the_cd, n+1, TRANSLOCATION);
                            rev.same_eds = TRUE;
                        }
                        else if ((n -= the_cd->n00_seds_trans) < the_cd->n00_deds1_trans){
                            rev = get_spec_00deds1_rev(the_cd, n+1, TRANSLOCATION);
                            rev.same_eds = FALSE;
                        }
                        else if((n -= the_cd->n00_deds1_trans) < the_cd->n00_deds2_trans){
                            rev = get_spec_00deds2_rev(the_cd, n+1, TRANSLOCATION);
                            rev.same_eds = FALSE;
                        }
                        else printf("in delta_ci_cg_reverse. d = 0. Problem.\n");   
                            //      printf("dd=0, translocation. \n"); 
                        *prob = p_dd0*p_trans/(double)the_cd->n_dd_trans[1];
                    }
                     

                }
                else{ // each inversion/translocation equally probably
                    n = (int)(drand(&rng_long)*n_dd0);
                    assert(n_dd0 == (the_cd->n00_seds + the_cd->n00_deds1 + the_cd->n00_deds2));
                    if(n < the_cd->n00_seds){
                        rev = get_spec_00seds_rev(the_cd, n+1, EITHER);
                        rev.same_eds = TRUE;
                    }
                    else if ((n -= the_cd->n00_seds) < the_cd->n00_deds1){
                        rev = get_spec_00deds1_rev(the_cd, n+1, EITHER);
                        rev.same_eds = FALSE;
                    }
                    else if((n -= the_cd->n00_deds1) < the_cd->n00_deds2){
                        rev = get_spec_00deds2_rev(the_cd, n+1, EITHER);
                        rev.same_eds = FALSE;
                    }
                    else printf("in delta_ci_cg_reverse. d = 0. Problem.\n");   
         
                    *prob = p_dd0/(double)n_dd0;
                }
            }
            rev.delta_c_i = 0; rev.delta_c_g = 0;
        }
        else{  // delta d = +1 step
            if(PREFER_INVERSIONS == TRUE){ // prefer inversions
                    //   double denom;
                if(USE_IPF_LAMBDA_RATIO){
                    p_inv = pow(lambda_ratio, IPF_EXPONENT)*the_cd->n_dd_inv[2];
                }
                else{p_inv = pow(1.0, IPF_EXPONENT)*the_cd->n_dd_inv[2]; }
                    //   printf("IPF^IPF_exponent: %9g ", pow(INVERSION_PREFER_FACTOR, IPF_EXPONENT));
                p_trans = 1.0*the_cd->n_dd_trans[2];
                denom = p_inv + p_trans;
                p_inv /= denom; p_trans /= denom; // normalize
                if(drand(&rng_long) < p_inv){ // choose from among the inversions
                    n = (int)(drand(&rng_long)*the_cd->n_dd_inv[2]);
                    if(n < the_cd->n_dcg_inv[2]){ // do a delta c_g = 1 step
                        rev = get_spec_dcg1_rev(the_cd, n+1, INVERSION);
                        rev.delta_c_i = 0; rev.delta_c_g = 1;
                    }
                    else{ // do a delta c_i = -1 step
                        n -= the_cd->n_dcg_inv[2];
                        rev = get_spec_dcim1_rev(the_cd, n+1, INVERSION);
                        rev.delta_c_i = -1; rev.delta_c_g = 0;
                    }
                    rev.same_eds = FALSE;
                        //    printf("p_dd1, n_dcg1, n_dcim1: %g  %i %i \n", p_dd1, n_dcg1, n_dcim1);
                    *prob = p_dd1*p_inv/(double)the_cd->n_dd_inv[2];
                        //   printf("dd=1, inversion. \n"); 
                }
                else { // choose from among the translocations
                    n = (int)(drand(&rng_long)*the_cd->n_dd_trans[2]);
                    if(n < the_cd->n_dcg_trans[2]){ // do a delta c_g = 1 step
                        rev = get_spec_dcg1_rev(the_cd, n+1, TRANSLOCATION);
                        rev.delta_c_i = 0; rev.delta_c_g = 1;
                    }
                    else{ // do a delta c_i = -1 step
                        n -= the_cd->n_dcg_trans[2];
                        rev = get_spec_dcim1_rev(the_cd, n+1, TRANSLOCATION);
                        rev.delta_c_i = -1; rev.delta_c_g = 0;
                    }
                    rev.same_eds = FALSE;
                        //    printf("p_dd1, n_dcg1, n_dcim1: %g  %i %i \n", p_dd1, n_dcg1, n_dcim1);
                    *prob = p_dd1*p_trans/(double)the_cd->n_dd_trans[2];
                        //   printf("dd=1, translocation. \n"); 
                } 
            }
            else{
                n = (int)(drand(&rng_long)*(n_dcg1 + n_dcim1));
                if(n < n_dcg1){ // do a delta c_g = 1 step
                    rev = get_spec_dcg1_rev(the_cd, n+1, EITHER);
                    rev.delta_c_i = 0; rev.delta_c_g = 1;
                }
                else{ // do a delta c_i = -1 step
                    n -= n_dcg1;
                    rev = get_spec_dcim1_rev(the_cd, n+1, EITHER);
                    rev.delta_c_i = -1; rev.delta_c_g = 0;
                }
                rev.same_eds = FALSE;
                    //    printf("p_dd1, n_dcg1, n_dcim1: %g  %i %i \n", p_dd1, n_dcg1, n_dcim1);
                *prob = p_dd1/(double)(n_dcg1 + n_dcim1);
            }
        }
    }
    if(PRINTSTUFF) printf(" n00s %4i %4i %4i %4i %4i; ", the_cd->n00_seds,  the_cd->n00_deds1, the_cd->n00_deds2,
                          distance - (n_dci1+n_dcgm1), n_dci1+n_dcgm1-the_cd->n00_seds);
    if(PRINTSTUFF) printf("%4i %4i\n", rev.delta_c_i, rev.delta_c_g);
    
    if(step_fcn(n_dcgm1 + n_dci1) == 0){
             //   printf("no delta_d = -1 steps available \n");
        rev.stuck = TRUE;
    }
    else rev.stuck = FALSE;
         // printf("rev.stuck: %i \n", rev.stuck);
  /*   printf("in delta_ci_cg_reverse_new. *prob: %g \n", *prob); */
/*     printf("rev: %i %i %i %i %i %i \n", rev.left, rev.right, rev.is_inversion, rev.same_eds, rev.delta_c_i, rev.delta_c_g); */
    return rev;
} // end delta_ci_cg_reverse_new

double
p_dnsw(int dnsw_index, int n)
{
        // returns the (rel) prob of doing a dnsw = dnsw_index-2 inv/trans
        // if there are n to choose from
    static double mins[5] = {10.0, 10.0, 10.0, 0.0, 0.0};
    static double slopes[5] = {1.0, 0.5, 0.25, 0.08, 0.04};
    static double maxes[5] = {1.0e20, 1.0e20, 1.0e20, 5.0, 2.5};
    double min = mins[dnsw_index];
    double lin = slopes[dnsw_index]*(double)n;
    double max = maxes[dnsw_index];
    double result;

    if(n == 0) return 0.0; // if there are no inv/trans of this dnsw, prob is zero!
    result = (lin > min)? lin: min;
    result = (result > max)? max: result;
        //  printf("dnsw_index, n, min, lin, max, result: %8i %8i  %12g %12g %12g  %12g \n", dnsw_index, n, min, lin, max, result);
    return result;
}

int
n_edge_pairs(Permutation* perm)
{
        // returns N+M choose 2 - i.e. the number of edge pairs
    int n = perm->n_mark + perm->n_chrom;

    return n*(n-1)/2;
}



int get_n_dci1_eds(int N_chrom, Cycle_element** start, int* n_dci0, int* n_dci1_inv, int* n_dci1_trans)
{
#define ALLCHROMOSOMES FALSE    
        // returns the number of delta c_internal = 1 inversions and translocations
        // on the end_delimited_segment starting at cycle element pointed to by start
        // and  sets start to ptr to beginning of next eds
    
   
    int n_dci1 = 0, n_temp1, n_temp2;
    Cycle_element* elem = *start;
        // int the_eds = elem->eds_in_cycle;
    int the_eds = elem->eds;

        // if only one edge on the eds return 0
    if((elem->next == *start) || (elem->next->eds != the_eds)) {
        *start = elem->next;
        *n_dci1_inv = 0; *n_dci1_trans = 0;
        return 0;
    }
    else{
        int* nr = (int*)chcalloc((size_t)N_chrom, sizeof(int));
        int* nl = (int*)chcalloc((size_t)N_chrom, sizeof(int));
        int* chroms;
        int n_chroms_visited = 0;
        int nr_tot = 0, nl_tot = 0;             
            //   int the_cycle = (*start)->cycle_number; //
        int M, i;
        int n_eds_elems = 0; // counts number of elems on eds
        int n_dci0_eds;
        int sum_r = 0, sum_l = 0, sum_rl = 0, sum_rr = 0, sum_ll = 0;
        int n_inv, n_trans_flip, n_trans_noflip;
       
        if(!ALLCHROMOSOMES) chroms = (int*)chcalloc((size_t)N_chrom, sizeof(int)); // holds chromosomes visited by eds
        do { // loop over elements in eds
            M = elem->chromosome;

            if(!ALLCHROMOSOMES){
                if((nl[M] == 0) && (nr[M] == 0)){
                    chroms[n_chroms_visited] = M;
                    n_chroms_visited++;
                }
            }
             
            if(elem->direction == RIGHT){ nr[M]++; nr_tot++;}
            else{ nl[M]++; nl_tot++; } // LEFT
            elem = elem->next; // go to next edge in cycle
            n_eds_elems++;

        } while((elem->eds == the_eds) && (elem != *start));
        *start = elem;
            // if !ALLCHROMOSOMES just consider chromosomes visited by the eds
            // flipped translocations
        n_temp1 = n_temp2 = 0;
        if(ALLCHROMOSOMES){
            for(i=0; i<N_chrom; i++){
                if(nr[i]+nl[i] > 0){
                        /*   n_dci1 += nr[i]*(nr_tot - nr[i]);  */
/*                     n_dci1 += nl[i]*(nl_tot - nl[i]); */
                    n_temp1 += nr[M]*nr[M] + nl[M]*nl[M];
                    n_temp2 += nr[M]*nl[M];
                }
            }
        }
        else{
            for(i=0; i<n_chroms_visited; i++){
                M = chroms[i];
                    //   n_dci1 += nr[M]*(nr_tot - nr[M]); 
                    //    n_dci1 += nl[M]*(nl_tot - nl[M]);
                n_temp1 += nr[M]*nr[M] + nl[M]*nl[M];
                n_temp2 += nr[M]*nl[M];
                sum_r += nr[M]; sum_l += nl[M];
                sum_rr += nr[M]*nr[M]; sum_ll += nl[M]*nl[M];
                sum_rl += nr[M]*nl[M];
            }
        }
        n_inv = sum_rl;
        n_trans_noflip = sum_r*sum_l - sum_rl;
        n_trans_flip = (sum_r*sum_r + sum_l*sum_l - sum_rr - sum_ll)/2;
       
        n_dci1 = nr_tot*nr_tot + nl_tot*nl_tot - n_temp1;
        assert(n_dci1 %2 == 0);
        n_dci1 /= 2; // loop double counted
        n_dci1 += nr_tot*nl_tot; // add inversions and unflipped translocations
        assert(n_dci1 == ((nr_tot + nl_tot)*(nr_tot + nl_tot) - n_temp1)/2);
        
        n_dci0_eds = ((nr_tot + nl_tot)*(nr_tot + nl_tot - 1))/2 - n_temp2;
        *n_dci0 += n_dci0_eds;

            //  printf("%i %i %i %i %i \n", n_inv, n_trans_noflip, n_trans_flip, n_inv + n_trans_noflip + n_trans_flip, n_dci1);
        free(nr); free(nl);
        if(!ALLCHROMOSOMES) free(chroms);
     
        *n_dci1_inv = n_inv;
        *n_dci1_trans = n_trans_flip + n_trans_noflip;
        assert(n_dci1 == n_inv+n_trans_flip+n_trans_noflip);
        return n_dci1;
    }
} // end of get_n_dci1_eds


Reversal
get_spec_dci1_rev(Cycle_decomposition* the_cd, int spec_rev, int type)
{
        // return the spec_rev'th Reversal with delta c_i = +1
        // spec_rev is 1,2,3...
    int i, n_dci1_type_old, the_eds;
    Cycle_element* start, * start_old, * edge1, * edge2, * next1;
    Cycle* eds;
    int n_dci0_seds;
    int n_dci1 = 0, n_dci1_inv = 0, n_dci1_trans = 0;
    int  n_dci1_inv_eds, n_dci1_trans_eds;
    int* n_dci1_type = NULL;

    if(type == EITHER){ n_dci1_type = &n_dci1; }
    else if(type == INVERSION){ n_dci1_type = &n_dci1_inv; }
    else if(type == TRANSLOCATION){ n_dci1_type = &n_dci1_trans; }
    else{ printf("in get_spec_dci1_rev. unknown type. \n");} 
  
    for(i=0, eds = the_cd->edss; i<the_cd->n_edss; i++, eds++){ // loop over edss
        start = eds->start_elem_ptr;
            
            //   for(j=0; j<n_eds; j++){ // loop over eds's in cycle
              
        start_old = start;
        n_dci1_type_old =  *n_dci1_type;
        n_dci1 += get_n_dci1_eds(the_cd->n_chrom, &start, &n_dci0_seds, &n_dci1_inv_eds, &n_dci1_trans_eds);
        n_dci1_inv += n_dci1_inv_eds;
        n_dci1_trans += n_dci1_trans_eds;
            //  printf("ndci1's : %i %i %i \n", n_dci1, n_dci1_inv, n_dci1_trans);
            //   n_dci1 += get_n_dci1_eds(the_cd->n_chrom, &start, &n_dci0_seds, type);
        if(spec_rev <= *n_dci1_type){ // go back to last eds done and find the specific reversal requested
            *n_dci1_type = n_dci1_type_old; // just set the one of interest (inv, trans, or either) back to beginning of eds
                //  n_dci1 = n_dci1_old;
                //   printf("type, *ndci1_type: %i %i \n", type, *n_dci1_type);
            the_eds = start_old->eds; //  printf("j, the_eds, n_eds: %i %i %i \n", j, the_eds, n_eds);
            edge1 = start_old;
                  
            for(edge1 = start_old; (edge1->next->eds == the_eds) && (edge1->next != start_old); edge1 = edge1->next){
                next1 = edge1->next;
                for(edge2 = next1; (edge2->eds == the_eds) && (edge2 != start_old); edge2 = edge2->next){
                    if(edge1->direction != edge2->direction){ // inversion or unflipped
                        n_dci1++;
                        if(edge1->chromosome != edge2->chromosome){ n_dci1_trans++;}
                        else { n_dci1_inv++; }
                        if(*n_dci1_type == spec_rev){
                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                        }
                    }
                    else { // same direction - flipped translocation if on distinct chromosomes
                        if(edge1->chromosome != edge2->chromosome){ // flipped translocation
                            n_dci1++; n_dci1_trans++;
                            if(*n_dci1_type == spec_rev){
                                return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                            }  
                        }
                    }
                }
            }
        }
    }
    printf("n_dci1, n_dci1_inv, n_dci1_trans: %i  %i  %i\n", n_dci1, n_dci1_inv, n_dci1_trans);
    printf("in get_spec_dci1_rev, couldn't find the specified rev\n");
} // end of get_spec_dci1_rev




Reversal
get_spec_00seds_rev(Cycle_decomposition* the_cd, int n_spec, int type)
{
   // returns the number of delta c_i = 1 , (and delta_c_g = 0) inversion/translocations
        // uses edss, not cycles
        // also counts number of delta c_i = 0, delta c_g = 0 inversion/translocations with both edges in same eds
    Cycle* eds;
    Cycle_element* edge1, * edge2;
    int i;
    int n_dci1 = 0;
    int n_00seds = 0, n_00seds_inv = 0, n_00seds_trans = 0; // use this to count number of delta c_i - delta_c_g = 0 inv/translocations on same edss
    int* n_00seds_type = NULL;
    
    if(type == EITHER){ n_00seds_type = &n_00seds; }
    else if(type == INVERSION){ n_00seds_type = &n_00seds_inv; }
    else if(type == TRANSLOCATION){ n_00seds_type = &n_00seds_trans; }
    else{ printf("in get_spec_00seds_rev. unknown type. \n");}
      
    for(i=0, eds = the_cd->edss; i<the_cd->n_edss; i++, eds++){
        if(eds->n_black > 1){
            edge1 = eds->start_elem_ptr;   
            do{
                for(edge2 = edge1->next; ((edge2->eds == i) && (edge2 != eds->start_elem_ptr)); edge2 = edge2->next){
                    if(edge1->chromosome == edge2->chromosome){
                        if(edge1->direction != edge2->direction) n_dci1++;
                        else {n_00seds++; n_00seds_inv++;}
                        if(*n_00seds_type == n_spec)
                           return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1); 
                    }
                    else { // translocation
                        n_dci1++; n_00seds++; n_00seds_trans++;
                    if(*n_00seds_type == n_spec)
                       return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, (edge1->direction != edge2->direction), the_cd->perm1);  
                        }
                }
                edge1 = edge1->next;
            }while((edge1->eds == i) && (edge1 != eds->start_elem_ptr));
        }
    }
    printf("in get_spec_00seds_rev. REached end without finding specified rev. \n");
}



Reversal
get_spec_00deds1_rev(Cycle_decomposition* the_cd, int n_spec, int type)
{
    // 
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j;
    int n_00deds1 = 0, n_00deds1_inv = 0, n_00deds1_trans = 0;
    int* n_00deds1_type = NULL; 

    if(type == EITHER){ n_00deds1_type = &n_00deds1; }
    else if(type == INVERSION){ n_00deds1_type = &n_00deds1_inv; }
    else if(type == TRANSLOCATION){ n_00deds1_type = &n_00deds1_trans; }
    else{ printf("in get_spec_00seds_rev. unknown type. \n");}
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle && (eds1->first_black != eds1->last_black)){
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle  && (eds2->first_black != eds2->last_black)){ 
                    if(eds1->first_black != eds2->first_black){ // bg + gb, opposite (same) directions required (unflipped)  for delta c_g = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        do{
                            edge2 = eds2->start_elem_ptr;
                            do{
                                if(edge1->chromosome == edge2->chromosome){ //inversion
                                    if(edge1->direction == edge2->direction){
                                        n_00deds1++; n_00deds1_inv++;
                                        if(*n_00deds1_type == n_spec){
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                    }
                                }
                                else {
                                    n_00deds1 += 1; n_00deds1_trans++; // exactly 1 of flipped, unflipped will be delta c_g = 1
                                    if(*n_00deds1_type == n_spec){
                                        if(edge1->direction == edge2->direction){
                                           return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                        else{
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                                        }
                                    }
                                }
                                edge2 = edge2->next;
                            }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                            edge1 = edge1->next;
                        }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                    }
                    else { // bg + bg or gb + gb, same (opp) directions required (unflipped) for delta c_g = 1 (0)
                        edge1 = eds1->start_elem_ptr;
                        do{
                            edge2 = eds2->start_elem_ptr;
                            do{
                               if(edge1->chromosome == edge2->chromosome){ //inversion
                                    if(edge1->direction != edge2->direction){
                                        n_00deds1++; n_00deds1_inv++;
                                        if(*n_00deds1_type == n_spec){
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                    }
                                }
                                else {
                                    n_00deds1 += 1; n_00deds1_trans += 1; // exactly 1 of flipped, unflipped will be delta c_g = 1
                                    if(*n_00deds1_type == n_spec){
                                        if(edge1->direction != edge2->direction){
                                           return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                        }
                                        else{
                                            return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                                        }
                                    }
                                }
                                edge2 = edge2->next;
                            }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                            edge1 = edge1->next;
                        }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                    }
                }
            }
        }
    }
    printf("In get_spec_00deds1_rev. Reached end without finding specified reversal.\n");
} // end of get_spec_00deds1_rev


Reversal
get_spec_00deds2_rev(Cycle_decomposition* the_cd, int n_spec, int type)
{
        // returns the number of delta c_g = -1 inversion/translocations
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j, nb1;
    int n_00deds2 = 0,  n_00deds2_inv = 0,  n_00deds2_trans = 0;
    int* n_00deds2_type  = NULL;

     if(type == EITHER){ n_00deds2_type = &n_00deds2; }
    else if(type == INVERSION){ n_00deds2_type = &n_00deds2_inv; }
    else if(type == TRANSLOCATION){ n_00deds2_type = &n_00deds2_trans; }
    else{ printf("in get_spec_00seds_rev. unknown type. \n");} 
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle){
            nb1 = eds1->first_black + eds1->last_black;
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle && ((nb1 + eds2->first_black + eds2->last_black) != 2)){ 
                    edge1 = eds1->start_elem_ptr;
                    do{
                        edge2 = eds2->start_elem_ptr;
                        do{
                            if(edge1->chromosome == edge2->chromosome){
                                n_00deds2++; n_00deds2_inv++;
                                if(*n_00deds2_type == n_spec)
                                    return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                            }
                            else { // translocation
                                n_00deds2++; n_00deds2_trans++;
                                if(*n_00deds2_type == n_spec)
                                    return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                n_00deds2++; n_00deds2_trans++;
                                if(*n_00deds2_type == n_spec)
                                    return make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1); 
                            }
                            edge2 = edge2->next;
                        }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                        edge1 = edge1->next;
                    }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                }
            }
        }
    }
    printf("In get_spec_00deds2_rev. Reached end without finding specified reversal.\n");
} // end of get_spec_00deds2_rev


int
get_n_dcgm1(Cycle_decomposition* the_cd, const Permutation* p2)
{
        // returns the number of delta c_g = -1 inversion/translocations
     Cycle* eds1, * eds2;
     Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcgm1 = 0;
    int n_dcgm1_inv = 0, n_dcgm1_trans = 0;
    int dnsw_inv[5] = {0,0,0,0,0}, dnsw_trans[5] = {0,0,0,0,0};
    int delta_n_switches_unflipped, delta_n_switches_flipped;
    int nswbef, nswaft, nswaftflip;
     int ncom_abcd, ncom_acbd, ncom_adbc;
     const Permutation* p1 = the_cd->perm1;
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle && (eds1->first_black == eds1->last_black)){
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle  && (eds2->first_black == eds2->last_black) && (eds1->first_black != eds2->first_black)){
                    edge1 = eds1->start_elem_ptr;
                    do{
                        edge2 = eds2->start_elem_ptr;
                        do{
                            delta_n_switches_unflipped = delta_n_switches(edge1, edge2, p2, &delta_n_switches_flipped, &nswbef, &nswaft, &nswaftflip);
                            if(edge1->chromosome == edge2->chromosome){ // count 1 inversion
                                n_dcgm1++; n_dcgm1_inv++;
                                    //  delta_n_switches_unflipped = delta_n_switches(edge1, edge2, p2, &delta_n_switches_flipped, &nswbef, &nswaft, &nswaftflip);
                                dnsw_inv[2 + delta_n_switches_unflipped]++;
                                    // if(nswbef == 2) the_cd->n_ddm1_nsw2before[2 + delta_n_switches_unflipped]++;
                            } 
                            else { // count 2 translocations
                                n_dcgm1 += 2; n_dcgm1_trans += 2;
                                dnsw_trans[2 + delta_n_switches_unflipped]++;
                                dnsw_trans[2 + delta_n_switches_flipped]++;
                                if(nswbef == 2){
                                      ncom_abcd = chroms_in_common(edge1, edge2, p1, p2, &ncom_acbd, &ncom_adbc);
                                          //      printf("in get_n_dcgm1. both trans. ncoms: orig, unfl, flip: %i %i %i \n", ncom_abcd, ncom_acbd, ncom_adbc); 
                                    the_cd->n_ddm1_nsw2before[2 + delta_n_switches_unflipped]++;
                                    the_cd->n_ddm1_nsw2before[2 + delta_n_switches_flipped]++;
                                    
                                   
                                    if(ncom_acbd < ncom_abcd) {
                                        the_cd->n_ddm1_ncomdecrease[2 + delta_n_switches_unflipped]++;
                                        if(ncom_acbd == 0) the_cd->n_ddm1_ncomdectozero[2 + delta_n_switches_unflipped]++;
                                    }
                                    else if(ncom_acbd == ncom_abcd) {
                                         the_cd->n_ddm1_ncomnochange[2 + delta_n_switches_unflipped]++;
                                    }
                                   
                                    if(ncom_adbc < ncom_abcd) {
                                        the_cd->n_ddm1_ncomdecrease[2 + delta_n_switches_flipped]++;
                                        if(ncom_adbc == 0) the_cd->n_ddm1_ncomdectozero[2 + delta_n_switches_flipped]++;
                                    }
                                     else if(ncom_adbc == ncom_abcd) {
                                         the_cd->n_ddm1_ncomnochange[2 + delta_n_switches_flipped]++;
                                    }
                                    
                                  /*   if(delta_n_switches_unflipped > 0 || delta_n_switches_flipped > 0){ */
/*                                         printf("nswbef, delta_n_switches, un, fflip: %i %i %i \n", nswbef, delta_n_switches_unflipped,  delta_n_switches_flipped); */
/*                                     } */
                                }
                            }
                            edge2 = edge2->next;
                        }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                        edge1 = edge1->next;
                    }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                }
            }
        }
    }
    the_cd->n_dcg_inv[0] = n_dcgm1_inv;
    the_cd->n_dcg_trans[0] = n_dcgm1_trans;

     for(j=0; j<5; j++){
        the_cd->n_dsw_inv_dcgm1[j] = dnsw_inv[j];
        the_cd->n_dsw_trans_dcgm1[j] = dnsw_trans[j];
        the_cd->n_dsw_inv_ddm1[j] += dnsw_inv[j];
        the_cd->n_dsw_trans_ddm1[j] += dnsw_trans[j];
        the_cd->n_dsw_dcgm1[j] = dnsw_inv[j] + dnsw_trans[j];
     }
    return n_dcgm1;
}


int
get_n_dci1(Cycle_decomposition* the_cd, const Permutation* p2)
{
        // returns the number of delta c_i = 1 , (and delta_c_g = 0) inversion/translocations
        // uses edss, not cycles
        // also counts number of delta c_i = 0, delta c_g = 0 inversion/translocations with both edges in same eds
    Cycle* eds;
    Cycle_element* edge1, * edge2;
    int i;
    int n_dci1 = 0;
    int n_dci1_inv = 0, n_dci1_trans = 0;
    int n00_seds_inv = 0, n00_seds_trans = 0;
    int n = 0; // use this to count number of delta c_i = delta_c_g = 0 inv/translocations on same edss
    int dnsw_inv[5] = {0,0,0,0,0}, dnsw_trans[5] = {0,0,0,0,0};
    int delta_n_switches_unflipped, delta_n_switches_flipped;
    int nswbef, nswaft, nswaftflip;
    const Permutation* p1 = the_cd->perm1;
    int ncom_abcd, ncom_acbd, ncom_adbc;
    
    
    for(i=0, eds = the_cd->edss; i<the_cd->n_edss; i++, eds++){
        if(eds->n_black > 1){
            edge1 = eds->start_elem_ptr;   
            do{
                for(edge2 = edge1->next; ((edge2->eds == i) && (edge2 != eds->start_elem_ptr)); edge2 = edge2->next){
                    if(edge1->chromosome == edge2->chromosome){ // inversion
                        if(edge1->direction != edge2->direction) {
                            n_dci1++; n_dci1_inv++; 
                            dnsw_inv[2 + delta_n_switches(edge1, edge2, p2, &delta_n_switches_flipped, &nswbef, &nswaft, &nswaftflip)]++;      
                        }
                        else {n++; n00_seds_inv++;}
                    }
                    else { // translocation
                        n_dci1++; n_dci1_trans++; n++; n00_seds_trans++;
                        delta_n_switches_unflipped = delta_n_switches(edge1, edge2, p2, &delta_n_switches_flipped, &nswbef, &nswaft, &nswaftflip);
                        if(edge1->direction != edge2->direction){ // unflipped is the one
                            dnsw_trans[2 + delta_n_switches_unflipped]++;
                            if(nswbef == 2) {
                                ncom_abcd = chroms_in_common(edge1, edge2, p1, p2, &ncom_acbd, &ncom_adbc);

                              
                                if(ncom_acbd < ncom_abcd) {
                                    the_cd->n_ddm1_ncomdecrease[2 + delta_n_switches_unflipped]++;
                                    if(ncom_acbd == 0) the_cd->n_ddm1_ncomdectozero[2 + delta_n_switches_unflipped]++;
                                }
                                 else if(ncom_acbd == ncom_abcd) {
                                         the_cd->n_ddm1_ncomnochange[2 + delta_n_switches_unflipped]++;
                                    }
                                  
                                    //   printf("in get_n_dci1.unflip trans: ncoms: orig, unfl, flip: %i %i %i \n", ncom_abcd, ncom_acbd, ncom_adbc);
                                the_cd->n_ddm1_nsw2before[2 + delta_n_switches_unflipped]++;
                              /*   if(delta_n_switches_unflipped > 0 || delta_n_switches_flipped > 0){ */
/*                                     printf("nswbef, delta_n_switches, un, fflip: %i %i %i \n", nswbef, delta_n_switches_unflipped,  delta_n_switches_flipped); */
/*                                 } */
                            }
                        }
                        else{ // flipped one is the one
                            dnsw_trans[2 + delta_n_switches_flipped]++;
                            if(nswbef == 2) {
                                ncom_abcd = chroms_in_common(edge1, edge2, p1, p2, &ncom_acbd, &ncom_adbc);

                              
                                if(ncom_adbc < ncom_abcd){
                                    the_cd->n_ddm1_ncomdecrease[2 + delta_n_switches_flipped]++;
                                    if(ncom_adbc == 0) the_cd->n_ddm1_ncomdectozero[2 + delta_n_switches_flipped]++;
                                }
                                 else if(ncom_acbd == ncom_adbc) {
                                         the_cd->n_ddm1_ncomnochange[2 + delta_n_switches_flipped]++;
                                    }
                                
                                    //  printf("in get_n_dci1. flippedtrans. ncoms: orig, unfl, flip: %i %i %i \n", ncom_abcd, ncom_acbd, ncom_adbc);
                                the_cd->n_ddm1_nsw2before[2 + delta_n_switches_flipped]++;
                              /*   if(delta_n_switches_unflipped > 0 || delta_n_switches_flipped > 0){ */
/*                                     printf("nswbef, delta_n_switches, un, flip: %i %i %i \n", nswbef, delta_n_switches_unflipped,  delta_n_switches_flipped); */
/*                                 } */ 
                            }
                        }
                            
                    } // one of each, delta c_i=0 and 00seds
                }
                edge1 = edge1->next;
            }while((edge1->eds == i) && (edge1 != eds->start_elem_ptr));
        }
    }
    
    the_cd->n00_seds = n;
    the_cd->n_dci_inv[2] = n_dci1_inv;
    the_cd->n_dci_trans[2] = n_dci1_trans;
    the_cd->n00_seds_inv = n00_seds_inv;
    the_cd->n00_seds_trans = n00_seds_trans;

    for(i=0; i<5; i++){
        the_cd->n_dsw_inv_dci1[i] = dnsw_inv[i];
        the_cd->n_dsw_trans_dci1[i] = dnsw_trans[i];
        the_cd->n_dsw_inv_ddm1[i] += dnsw_inv[i];
        the_cd->n_dsw_trans_ddm1[i] += dnsw_trans[i];
        the_cd->n_dsw_dci1[i] = dnsw_inv[i] +  dnsw_trans[i];
    }
    return n_dci1;
} // end of get_n_dci1
        

Reversal
get_spec_dci1_rev_new(Cycle_decomposition* the_cd, int spec_rev, int type, int spec_dsw_index)
{
        // return the spec_rev'th Reversal with delta c_i = +1
        // spec_rev is 1,2,3...
    int i, n_dci1_type_old, the_eds;
    Cycle_element* start, * start_old, * edge1, * edge2, * next1;
    Cycle* eds;
    int n_dci0_seds;
    int n_dci1 = 0, n_dci1_inv = 0, n_dci1_trans = 0;
    int  n_dci1_inv_eds, n_dci1_trans_eds;
    int* n_dci1_type  = NULL;
    int n_dci1_dsw[5] = {0};
//    int n_dci1_dsw_eds[5] = {0};
    int dnsw, dnsw_flipped;
    Reversal the_rev;
        //    int dsw1, dsw1_flipped, dsw2, dsw2_flipped;
      int nswbef, nswaft, nswaftflip; 

    if(type == EITHER){ n_dci1_type = &n_dci1; }
    else if(type == INVERSION){ n_dci1_type = &n_dci1_inv; }
    else if(type == TRANSLOCATION){ n_dci1_type = &n_dci1_trans; }
    else{ printf("in get_spec_dci1_rev. unknown type. \n");} 
  
    for(i=0, eds = the_cd->edss; i<the_cd->n_edss; i++, eds++){ // loop over edss
        start = eds->start_elem_ptr;
            
            //   for(j=0; j<n_eds; j++){ // loop over eds's in cycle
              
        start_old = start;
        n_dci1_type_old =  *n_dci1_type;
        n_dci1 += get_n_dci1_eds(the_cd->n_chrom, &start, &n_dci0_seds, &n_dci1_inv_eds, &n_dci1_trans_eds);
        n_dci1_inv += n_dci1_inv_eds;
        n_dci1_trans += n_dci1_trans_eds;
            //    printf("ndci1's : %i %i %i \n", n_dci1, n_dci1_inv, n_dci1_trans);
            //   n_dci1 += get_n_dci1_eds(the_cd->n_chrom, &start, &n_dci0_seds, type);
        if(spec_rev <= *n_dci1_type){ // go back to last eds done and find the specific reversal requested
            *n_dci1_type = n_dci1_type_old; // just set the one of interest (inv, trans, or either) back to beginning of eds
                //  n_dci1 = n_dci1_old;
                //   printf("type, *ndci1_type: %i %i \n", type, *n_dci1_type);
            the_eds = start_old->eds; //  printf("j, the_eds, n_eds: %i %i %i \n", j, the_eds, n_eds);
            edge1 = start_old;
                  
            for(edge1 = start_old; (edge1->next->eds == the_eds) && (edge1->next != start_old); edge1 = edge1->next){
                next1 = edge1->next;
             
                for(edge2 = next1; (edge2->eds == the_eds) && (edge2 != start_old); edge2 = edge2->next){
                        //    printf("in get_spec_dci1... before delta_n_sw... \n");
                    
                    dnsw = delta_n_switches(edge1, edge2, the_cd->perm2, &dnsw_flipped, &nswbef, &nswaft, &nswaftflip);
                
                    if(edge1->direction != edge2->direction){ // inversion or unflipped
                        n_dci1++;
                            //  dnsw = delta_n_switches(edge1, edge2, perm2, &dnsw_flipped);
                        n_dci1_dsw[dnsw+2]++;
                        if(edge1->chromosome != edge2->chromosome){ n_dci1_trans++;}
                        else { n_dci1_inv++; }
                        if(n_dci1_dsw[spec_dsw_index] == spec_rev){
                          /*   printf("inv or unflipped trans\n"); */
/*                             printf("in get_spec_rev_dci1_new. chroms a,b,c,d: %i %i %i %i    dnsw,noflip, flip: %i %i \n", */
/*                                    the_cd->perm2->p[edge1->n1].chromosome, */
/*                                    the_cd->perm2->p[edge1->n2].chromosome, */
/*                                    the_cd->perm2->p[edge2->n1].chromosome, */
/*                                    the_cd->perm2->p[edge2->n2].chromosome, */
/*                                    dnsw, dnsw_flipped); */
                            the_rev =  make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                          /*    dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                              dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                              if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped) */
/*                                         printf("in get_spec_dci1_rev_new. dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                                spec_dsw_index - 2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                                        
                                     /*   printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                        the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                 printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                        the_rev.numbs[2], the_rev.numbs[3]);  */  
                            return the_rev;
                        }
                    }
                    else { // same direction - flipped translocation if on distinct chromosomes
                        if(edge1->chromosome != edge2->chromosome){ // flipped translocation
                            n_dci1++; n_dci1_trans++;
                            n_dci1_dsw[dnsw_flipped+2]++;
                            if(n_dci1_dsw[spec_dsw_index] == spec_rev){
                              /*   printf("flipped trans\n"); */
/*                                 printf("in get_spec_rev_dci1_new. chroms a,b,c,d: %i %i %i %i    dnsw,noflip, flip: %i %i \n", */
/*                                        the_cd->perm2->p[edge1->n1].chromosome, */
/*                                        the_cd->perm2->p[edge1->n2].chromosome, */
/*                                        the_cd->perm2->p[edge2->n1].chromosome, */
/*                                        the_cd->perm2->p[edge2->n2].chromosome, */
/*                                        dnsw, dnsw_flipped); */
                                the_rev =   make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                            /*      dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                  dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                                    if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped) */
/*                                         printf("in get_spec_dci1_rev_new. dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                                spec_dsw_index - 2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                                      /*     printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                        the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                 printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                        the_rev.numbs[2], the_rev.numbs[3]); */
                                return the_rev;
                            }  
                        }
                    }
                    printf("n_dci1, total, inv, trans: %i %i %i \n", n_dci1, n_dci1_inv, n_dci1_trans);
                }
            }
        }
    }
    printf("n_dci1, n_dci1_inv, n_dci1_trans: %i  %i  %i\n", n_dci1, n_dci1_inv, n_dci1_trans);
    printf("n_dci1_dsw[spec_dsw_index] : %i  spec_rev: %i \n", n_dci1_dsw[spec_dsw_index], spec_rev);
    printf("in get_spec_dci1_rev, couldn't find the specified rev\n");
} // end of get_spec_dci1_rev_new

Reversal
get_spec_dci1_rev_new1(Cycle_decomposition* the_cd, int spec_n, int type, int spec_dnsw_index)
{
        // returns the number of delta c_i = 1 , (and delta_c_g = 0) inversion/translocations
        // uses edss, not cycles
        // also counts number of delta c_i = 0, delta c_g = 0 inversion/translocations with both edges in same eds
        // type is EITHER, INVERSION or TRANSLOCATION
    Cycle* eds;
    Cycle_element* edge1, * edge2;
    int i;
    int n_dci1 = 0;
    int n_dci1_inv = 0, n_dci1_trans = 0;
    int n00_seds_inv = 0, n00_seds_trans = 0;
    int n = 0; // use this to count number of delta c_i = delta_c_g = 0 inv/translocations on same edss
    int dnsw_inv[5] = {0,0,0,0,0}, dnsw_trans[5] = {0,0,0,0,0};
    int dnsw_both[5] = {0,0,0,0,0};
    int* dnsw_type = NULL ;
    int delta_n_switches_unflipped, delta_n_switches_flipped;
    Reversal the_rev;
        //   int dsw1, dsw1_flipped, dsw2, dsw2_flipped;
     int nswbef, nswaft, nswaftflip;
     
        // int dnsw;
    if(type == EITHER) { dnsw_type = dnsw_both; }
    else if(type == INVERSION) { dnsw_type = dnsw_inv;  }
    else if(type == TRANSLOCATION) { dnsw_type = dnsw_trans; }
    else printf("in get_spec_dci1_rev_new1. Unknown type\n"); 
        
    
    for(i=0, eds = the_cd->edss; i<the_cd->n_edss; i++, eds++){
        if(eds->n_black > 1){
            edge1 = eds->start_elem_ptr;   
            do{
                for(edge2 = edge1->next; ((edge2->eds == i) && (edge2 != eds->start_elem_ptr)); edge2 = edge2->next){
                    if(edge1->chromosome == edge2->chromosome){ // inversion
                        if(edge1->direction != edge2->direction) {
                            n_dci1++; n_dci1_inv++;
                            delta_n_switches_unflipped = delta_n_switches(edge1, edge2, the_cd->perm2, &delta_n_switches_flipped, &nswbef, &nswaft, &nswaftflip);
                            dnsw_inv[2 + delta_n_switches_unflipped]++; dnsw_both[2 + delta_n_switches_unflipped]++;
                            if(dnsw_type[spec_dnsw_index] == spec_n) {
                                the_rev = make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                             /*    dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                 dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
                                
/*                                 if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped)                     */
/*                                     printf("in get_spec_dci1_rev_new1. dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                            spec_dnsw_index-2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                                    /*   printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                        the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                 printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                        the_rev.numbs[2], the_rev.numbs[3]); */
                                return the_rev;
                            }
                        }
                        else {n++; n00_seds_inv++;}
                    }
                    else { // translocation
                        n_dci1++; n_dci1_trans++; n++; n00_seds_trans++;
                        delta_n_switches_unflipped = delta_n_switches(edge1, edge2, the_cd->perm2, &delta_n_switches_flipped, &nswbef, &nswaft, &nswaftflip);
                        if(edge1->direction != edge2->direction){ // unflipped is the one
                            dnsw_trans[2 + delta_n_switches_unflipped]++;  dnsw_both[2 + delta_n_switches_unflipped]++;
                            if(dnsw_type[spec_dnsw_index] == spec_n) {
                                the_rev = make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                              /*   dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                 dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                                   if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped)     */
/*                                 printf("in get_spec_dci1_rev_new1. dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                        spec_dnsw_index-2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                            /*     printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                        the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                 printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                        the_rev.numbs[2], the_rev.numbs[3]);   */
                                return the_rev;
                            }
                        }
                        else{ // flipped one is the one
                            dnsw_trans[2 + delta_n_switches_flipped]++; dnsw_both[2 + delta_n_switches_flipped]++;
                            if(dnsw_type[spec_dnsw_index] == spec_n) {
                                the_rev = make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                              /*   dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                 dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                                 if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped)       */
/*                                 printf("in get_spec_dci1_rev_new1. dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                        spec_dnsw_index-2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                             /*    printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                        the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                 printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                        the_rev.numbs[2], the_rev.numbs[3]);  */
                                return the_rev;
                            }
                        }       
                    }
                        //  printf("spec_dnsw_index: %i; dnsw_both/inv/trans[spec_dnsw_index]: %i %i %i; dnsw_type[spec...] %i, spec_n: %i \n",
                        //        spec_dnsw_index, dnsw_both[spec_dnsw_index], dnsw_inv[spec_dnsw_index], dnsw_trans[spec_dnsw_index], dnsw_type[spec_dnsw_index], spec_n);
                }
                edge1 = edge1->next;
            }while((edge1->eds == i) && (edge1 != eds->start_elem_ptr));
        }
    }
    for(i=0; i<5; i++){ printf("%i   ", the_cd->n_dsw_inv_dci1[i]); } printf("\n");
    for(i=0; i<5; i++){ printf("%i   ", the_cd->n_dsw_trans_dci1[i]); } printf("\n");
    printf("spec_dnsw_index: %i; dnsw_both/inv/trans[spec_dnsw_index]: %i %i %i; dnsw_type[spec...] %i, spec_n: %i \n",
                       spec_dnsw_index, dnsw_both[spec_dnsw_index], dnsw_inv[spec_dnsw_index], dnsw_trans[spec_dnsw_index], dnsw_type[spec_dnsw_index], spec_n);            
    printf("in get_spec_dci1_rev_new1.  Reached end without finding specified reversal.\n");
} // end of get_spec_dci1_rev_new1


Reversal
get_spec_dcgm1_rev_new(Cycle_decomposition* the_cd, int spec_n, int type, int spec_dsw_index)
{
        // returns the number of delta c_g = -1 inversion/translocations
    Cycle* eds1, * eds2;
    Cycle_element* edge1, * edge2;
    int i, j;
    int n_dcgm1 = 0,  n_dcgm1_inv = 0, n_dcgm1_trans = 0 ;
    int* n_dcgm1_dsw_type  = NULL;
    int n_dcgm1_dsw_both[5] = {0};
    int n_dcgm1_dsw_inv[5] = {0};
    int n_dcgm1_dsw_trans[5] = {0};
    int dnsw, dnsw_flipped;
    Reversal the_rev;
        //  int dsw1, dsw2, dsw1_flipped, dsw2_flipped;
    int nswbef, nswaft, nswaftflip;

    if(type == EITHER){ n_dcgm1_dsw_type = n_dcgm1_dsw_both; }
    else if(type == INVERSION){ n_dcgm1_dsw_type = n_dcgm1_dsw_inv; }
    else if(type == TRANSLOCATION){ n_dcgm1_dsw_type = n_dcgm1_dsw_trans; }
    else{ printf("in get_spec_dcgm1_rev_new. unknown type. \n");}
    
    for(i=0, eds1 = the_cd->edss; i<the_cd->n_edss; i++, eds1++){
        if(eds1->is_end_cycle && (eds1->first_black == eds1->last_black)){
            for(j = i+1, eds2 = the_cd->edss + j; j < the_cd->n_edss; j++, eds2++){
                if(eds2->is_end_cycle  && (eds2->first_black == eds2->last_black) && (eds1->first_black != eds2->first_black)){
                    edge1 = eds1->start_elem_ptr;
                    do{
                        edge2 = eds2->start_elem_ptr;
                        do{
                            if(edge1->chromosome == edge2->chromosome){ //inversion
                                n_dcgm1++; n_dcgm1_inv++;
                                dnsw = delta_n_switches(edge1, edge2, the_cd->perm2, &dnsw_flipped, &nswbef, &nswaft, &nswaftflip);
                                n_dcgm1_dsw_both[dnsw+2]++; n_dcgm1_dsw_inv[dnsw+2]++;
                                if(n_dcgm1_dsw_type[spec_dsw_index] == spec_n){ // inversion
                                  
                                  /*   printf("in get_spec_rev_dcgm1_new. Inversion. chroms a,b,c,d: %i %i %i %i    dnsw,noflip, flip: %i %i \n", */
/*                                            the_cd->perm2->p[edge1->n1].chromosome, */
/*                                            the_cd->perm2->p[edge1->n2].chromosome, */
/*                                            the_cd->perm2->p[edge2->n1].chromosome, */
/*                                            the_cd->perm2->p[edge2->n2].chromosome, */
/*                                            dnsw, dnsw_flipped); */
                                    
                                    the_rev = make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                  /*   dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                     dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                                    if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped)        */
/*                                     printf("in get_spec_dcgm1... dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                            spec_dsw_index - 2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                                 /*    printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                            the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                     printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                            the_rev.numbs[2], the_rev.numbs[3]); */
                                    return the_rev;
                                }
                            }
                            else { // translocation
                                n_dcgm1 += 2; n_dcgm1_trans += 2;
                                dnsw = delta_n_switches(edge1, edge2, the_cd->perm2, &dnsw_flipped, &nswbef, &nswaft, &nswaftflip);
                                n_dcgm1_dsw_both[dnsw+2]++; n_dcgm1_dsw_trans[dnsw+2]++;
                                    //    printf("dnsw, n_dcgm1_dsw_both[dnsw+2] : %i %i \n", dnsw, n_dcgm1_dsw_both[dnsw+2]);
                                if(n_dcgm1_dsw_type[spec_dsw_index] == spec_n){ // unflipped translocation
                                  /*   printf("in get_spec_rev_dcgm1_new. unfltrans. chroms a,b,c,d: %i %i %i %i    dnsw,noflip, flip: %i %i \n", */
/*                                            the_cd->perm2->p[edge1->n1].chromosome, */
/*                                            the_cd->perm2->p[edge1->n2].chromosome, */
/*                                            the_cd->perm2->p[edge2->n1].chromosome, */
/*                                            the_cd->perm2->p[edge2->n2].chromosome, */
/*                                            dnsw, dnsw_flipped); */
                                    the_rev =  make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, FALSE, the_cd->perm1);
                                 /*    dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                     dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                                      if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped)      */
/*                                     printf("in get_spec_dcgm1... dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                            spec_dsw_index - 2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                                /*     printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                            the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                     printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                            the_rev.numbs[2], the_rev.numbs[3]); */
                                    return the_rev;
                                }
                                else {
                                    n_dcgm1_dsw_both[dnsw_flipped+2]++; n_dcgm1_dsw_trans[dnsw_flipped+2]++;
                                    if (n_dcgm1_dsw_type[spec_dsw_index] == spec_n){ // flipped translocation
                                    
                                      /*   printf("in get_spec_rev_dcgm1_new. fliptrans. chroms a,b,c,d: %i %i %i %i    dnsw,noflip, flip: %i %i \n", */
/*                                                the_cd->perm2->p[edge1->n1].chromosome, */
/*                                                the_cd->perm2->p[edge1->n2].chromosome, */
/*                                                the_cd->perm2->p[edge2->n1].chromosome, */
/*                                                the_cd->perm2->p[edge2->n2].chromosome, */
/*                                                dnsw, dnsw_flipped); */
                                        the_rev =  make_rev(the_cd->n_mark + the_cd->n_chrom, edge1, edge2, TRUE, the_cd->perm1);
                                      /*   dsw1 = delta_n_switches1(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw1_flipped); */
/*                                         dsw2 = delta_n_switches2(the_rev.marker_numbers, the_cd->perm1, the_cd->perm2, &dsw2_flipped); */
/*                                          if(dsw1 != dsw2 || dsw1_flipped != dsw2_flipped)      */
/*                                         printf("in get_spec_dcgm1... dnsw's: %i   %i %i   %i %i flipped: %i\n", */
/*                                                spec_dsw_index - 2, dsw1, dsw1_flipped, dsw2, dsw2_flipped, the_rev.flip_chromosome); */
                                    /*     printf("%3i %3i %3i %3i \n", the_rev.marker_numbers[0], the_rev.marker_numbers[1], */
/*                                                the_rev.marker_numbers[2], the_rev.marker_numbers[3]); */
/*                                         printf("%3i %3i %3i %3i \n", the_rev.numbs[0], the_rev.numbs[1], */
/*                                                the_rev.numbs[2], the_rev.numbs[3]); */
                                        return the_rev;
                                    }
                                }
                            }
                            edge2 = edge2->next;
                        }while((edge2->eds == j) && (edge2 != eds2->start_elem_ptr));
                        edge1 = edge1->next;
                    }while((edge1->eds == i) && (edge1 != eds1->start_elem_ptr));
                }
            }
        }
    }
    for(i=0; i<5; i++){ printf("%i   ", the_cd->n_dsw_inv_dcgm1[i]); } printf("\n");
    for(i=0; i<5; i++){ printf("%i   ", the_cd->n_dsw_trans_dcgm1[i]); } printf("\n");
    printf("spec_dsw_index: %i; dnsw_both/inv/trans[spec_dsw_index]: %i %i %i; dsw_type[spec...] %i, spec_n: %i \n",
           spec_dsw_index,  n_dcgm1_dsw_both[spec_dsw_index],  n_dcgm1_dsw_inv[spec_dsw_index],
           n_dcgm1_dsw_trans[spec_dsw_index], n_dcgm1_dsw_type[spec_dsw_index], spec_n); 
    printf("in get_spec_dcgm1_rev_new. Reached end without finding specified reversal.\n");
} // end of get_spec_dcgm1_rev_new


int
chroms_in_common(Cycle_element* edge1, Cycle_element* edge2, const Permutation* p1, const Permutation* p2,
                 int* ncom_acbd, int* ncom_adbc)
{
        // given two edges, we want to know how many chromosomes (of genome p2) are represented on both chromosomes (i.e. the two to which the
        // two edges belong. This is ncom_abcd. WE also want to know the same thing for the rearranged genomes that would be produced by translocations
        // involvin the two edges. These are ncom_acbd (unflipped) and ncom_adbc (flipped).
    int achroms[MAX_N_CHROMOSOMES] = {0}, bchroms[MAX_N_CHROMOSOMES] = {0}, cchroms[MAX_N_CHROMOSOMES] = {0}, dchroms[MAX_N_CHROMOSOMES] = {0};
    int a,b,c,d;
    Marker_end* e;
    int ncom_abcd = 0;
    int i, n, chrom;
        //  Marker_end* a, * b, * c, * d;

    *ncom_acbd = 0; *ncom_adbc = 0;

    if(SHORTCIRCUITCHROMSINCOMMON) return 0;
    
    if(edge1->direction == LEFT){
        a = edge1->n1;
        b = edge1->n2;
    }
    else{
        a = edge1->n2;
        b = edge1->n1;
    }
    if(edge2->direction == LEFT){
        c = edge2->n1;
        d = edge2->n2;
    }
    else{
        c = edge2->n2;
        d = edge2->n1;
    }

    for(e=&(p1->p[a]); e->is_end_marker == FALSE; e = (e->other)->black){
        n = e->markend_numb;
        chrom = p2->p[n].chromosome; // chromosome of marker in p2
        achroms[chrom] = TRUE;
    }
    for(e=&(p1->p[b]); e->is_end_marker == FALSE; e = (e->other)->black){
        n = e->markend_numb;
        chrom = p2->p[n].chromosome; // chromosome of marker in p2
        bchroms[chrom] = TRUE;
    }
    for(e=&(p1->p[c]); e->is_end_marker == FALSE; e = (e->other)->black){
        n = e->markend_numb;
        chrom = p2->p[n].chromosome; // chromosome of marker in p2
        cchroms[chrom] = TRUE;
    }
    for(e=&(p1->p[d]); e->is_end_marker == FALSE; e = (e->other)->black){
        n = e->markend_numb;
        chrom = p2->p[n].chromosome; // chromosome of marker in p2
        dchroms[chrom] = TRUE;
    }

    for(i=0; i<MAX_N_CHROMOSOMES; i++){
        ncom_abcd += (achroms[i] || bchroms[i]) && (cchroms[i] || dchroms[i]) ;
        *ncom_acbd += (achroms[i] || cchroms[i]) && (bchroms[i] || dchroms[i]) ;
        *ncom_adbc += (achroms[i] || dchroms[i]) && (bchroms[i] || cchroms[i]) ;
    }
        //  printf("in chroms_in_common. ncoms: %i %i %i \n", ncom_abcd, *ncom_acbd, *ncom_adbc);
    return ncom_abcd;
}


void print_oxford_grid(FILE* fp, Permutation* perm1, Permutation* perm2)
{
    int n_chrom1 = perm1->n_chrom, n_chrom2 = perm2->n_chrom;
    int** ox = alloc_int2d(n_chrom1, n_chrom2);
    int i, j;
    Marker_end* m1, * m2;

    for(i=0; i<perm1->n_mark; i++){
        m1 = perm1->p + (2*i + 1);
        m2 = perm2->p + (2*i + 1);
        if(m1->is_end_marker){printf("int print_oxford_grid. m1 is end marker. \n"); getchar(); }
        if(m2->is_end_marker){printf("int print_oxford_grid. m2 is end marker. \n"); getchar(); }
        printf("m1, m2->markend_numb: %i %i \n", m1->markend_numb, m2->markend_numb);
        ox[m1->chromosome][m2->chromosome] += 1;
    }

    for(i=0; i<n_chrom1; i++){
        for(j=0; j<n_chrom2; j++){
            fprintf(fp, "%4i ", ox[i][j]); 
        } fprintf(fp, "\n");
    } fprintf(fp, "\n");

    free_int2d(n_chrom1, ox);
}

// end of perm.c
