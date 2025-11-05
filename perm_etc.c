// perm_etc.c

// stable stuff from perm.c Oct 31 02

#include "grit.h"
#include "params.h"
#include "structs.h"
#include "path_etc.h"
#include "perm_etc.h"
#include "exts.h"
#include "gr_misc.h"
#include "perm.h"

// get_n_inv etc.

int
get_n_inv(Permutation* perm)
{
        // returns the number of inversions, i.e.
        // sum over chromosomes of m choose 2 where m is number of edges on chromosome (i.e. # markers plus 1)
        //  int n_markends_on_chrom[MAXNCHROM] = {0};
  int* n_markends_on_chrom = (int*)chcalloc(perm->n_chrom, sizeof(int));
    int i, chrom, max_chrom = -1, nedges;
    int n_inv = 0; // count the number of inversions

    for(i=0; i<perm->n_chrom; i++){
        n_markends_on_chrom[i] = 0;
    }
    for(i=0; i<2*(perm->n_mark + perm->n_chrom); i++){
        chrom = perm->p[i].chromosome;
        n_markends_on_chrom[chrom]++;
        if(chrom > max_chrom) max_chrom = chrom;
    }

    assert(max_chrom <= perm->n_chrom);
    for(i=0; i<perm->n_chrom; i++){
        nedges = n_markends_on_chrom[i]; assert((nedges % 2) == 0);
        nedges /= 2;
        n_inv += nedges*(nedges-1)/2;
    }
    free(n_markends_on_chrom);
    return n_inv;
}


void
set_Ai_At(Permutation* perm)
{
    double sum_l_sqr = 0.0;
    double sum_chrom_l = 0.0, sum_chrom_l_sqr = 0.0;
    double sum_li_lj_inv = 0.0, sum_li_lj_trans = 0.0;
    double sum_ls = 0.0;
    double L_g = 0.0;
    int i; //, NpM = perm->n_mark + perm->n_chrom;
    Marker_end* e1, * e2;
    double l, area;
    int chrom;
    double* chrom_lengths = (double*)chcalloc(perm->n_chrom, sizeof(double));
    int stop = FALSE;

    /*  for(i=0; i<perm->n_chrom; i++){ */
    /* 	fprintf(stderr, "i, chrom_lengths[i]: %ld %g\n", i, chrom_lengths[i]); */
    /* } */
   
    for(e1 = perm->p; e1 != NULL; e1 = e2->other){
        e2 = e1->black;
        l = e2->distance - e1->distance;
            //     printf("UNKNOWN, e1->dist, e2->dist: %g %g %g \n", (double)UNKNOWN, e1->distance, e2->distance);
        if(e1->distance == UNKNOWN || e2->distance == UNKNOWN) {
	  fprintf(stderr, "In set_Ai_At. setting to UNKNOWN\n"); getchar();
            perm->A_i = UNKNOWN; perm->A_t = UNKNOWN; return;
        }
        if(l < 0.0 && DBDIST) {
            printf("in set_Ai_At. l negative : %g \n", l);
            stop = TRUE;
        }
        sum_l_sqr += l*l;
        L_g += l;
            //   printf("%g %g \n", first_l, sum_l_sqr);
        chrom = e1->chromosome;
        chrom_lengths[chrom] += l;
        sum_ls += l;
	// fprintf(stderr, "chrom, l, sum_ls, chrom_lengths[chrom]: %d %g    %g %g\n", chrom, l,  sum_ls, chrom_lengths[chrom]);
    }
    for(i=0; i<perm->n_chrom; i++){
        sum_chrom_l_sqr += chrom_lengths[i]*chrom_lengths[i];
        sum_chrom_l += chrom_lengths[i];
        perm->chromosome_lengths[i] = chrom_lengths[i];
	//fprintf(stderr, "i, chrom_lengths[i]: %ld %g\n", i, chrom_lengths[i]);
    }
    //printf("XCVB   %g %g %g %g\n", sum_ls, sum_chrom_l, sum_chrom_l_sqr, L_g);
    perm->A_t = L_g*L_g - sum_chrom_l_sqr;
    perm->A_i = 0.5*(sum_chrom_l_sqr - sum_l_sqr);
    if(FALSE || stop) {printf("in set_Ai_At. %g %g\n", perm->A_i, perm->A_t ); print_perm(stdout, perm); }// getchar();}
    
        // sum_l_sqr *= 0.5;
        //  printf("sum_l_sqr, sum_chrom_l_sqr, L_g^2: %g  %g  %g \n", sum_l_sqr, sum_chrom_l_sqr, L_g*L_g);
        //   printf("%g %g %g %g \n", sum_l_sqr, 2.0*perm->A_i, perm->A_t, L_g*L_g);
    free(chrom_lengths);
} // end of set_Ai_At


int
check_Ai_At(Permutation* perm)
{
    double sum_l_sqr = 0.0;
    double sum_chrom_l = 0.0, sum_chrom_l_sqr = 0.0;
    double sum_li_lj_inv = 0.0, sum_li_lj_trans = 0.0;
    double sum_ls = 0.0;
    double L_g = 0.0;
    int i; //, NpM = perm->n_mark + perm->n_chrom;
    Marker_end* e1, * e2;
    double l, area;
    int chrom;
    double* chrom_lengths = (double*)chcalloc(perm->n_chrom, sizeof(double));
    int stop = FALSE;
    int result = TRUE;
   
    for(e1 = perm->p; e1 != NULL; e1 = e2->other){
        e2 = e1->black;
        l = e2->distance - e1->distance;
            //     printf("UNKNOWN, e1->dist, e2->dist: %g %g %g \n", (double)UNKNOWN, e1->distance, e2->distance);
        if(e1->distance == UNKNOWN || e2->distance == UNKNOWN) {
            perm->A_i = UNKNOWN; perm->A_t = UNKNOWN; return result; // 
        }
        if(l < 0.0) {
            printf("in set_Ai_At. l negative : %g \n", l);
            stop = TRUE;
        }
        sum_l_sqr += l*l;
        L_g += l;
            //   printf("%g %g \n", first_l, sum_l_sqr);
        chrom = e1->chromosome;
        chrom_lengths[chrom] += l;
        sum_ls += l;
    }
    for(i=0; i<perm->n_chrom; i++){
        sum_chrom_l_sqr += chrom_lengths[i]*chrom_lengths[i];
        sum_chrom_l += chrom_lengths[i];
        if(perm->chromosome_lengths[i] != chrom_lengths[i]){
            printf("in check_Ai_At: i, perm->chromosome_lengths[i], check: %d %g %g \n", i, perm->chromosome_lengths[i], chrom_lengths[i]);
            result = FALSE;
        }    
        
    }
        //   printf("%g %g %g \n", sum_ls, sum_chrom_l, L_g);
    if(perm->A_t != L_g*L_g - sum_chrom_l_sqr) {
        printf("in check_Ai_At: At, Atcheck: %g %g \n", perm->A_t,  L_g*L_g - sum_chrom_l_sqr);
        result = FALSE;
    }
    if(perm->A_i != 0.5*(sum_chrom_l_sqr - sum_l_sqr)) {
        printf("in check_Ai_At: Ai, Aicheck: %g %g \n", perm->A_i,  0.5*(sum_chrom_l_sqr - sum_l_sqr));
        result = FALSE;
    }
    if(stop) {printf("in check_Ai_At. \n"); print_perm(stdout, perm); getchar();}
    
        // sum_l_sqr *= 0.5;
        //  printf("sum_l_sqr, sum_chrom_l_sqr, L_g^2: %g  %g  %g \n", sum_l_sqr, sum_chrom_l_sqr, L_g*L_g);
        //   printf("%g %g %g %g \n", sum_l_sqr, 2.0*perm->A_i, perm->A_t, L_g*L_g);
    free(chrom_lengths);
    return result;
}       // end of check_Ai_At


int
get_n_inv1(Permutation* perm)
{
  int markers_on_chrs[MAX_N_CHROMOSOMES] = {0}; // 
  int i, n_inv = 0;
  Marker_end* m = perm->p;

  m++;
  for(i=1; i<=perm->n_mark; i++, m += 2){
    if(m->is_end_marker){
      printf("in get_n_inv; marker is end marker. i= %i \n", i);
    }
    else{
      assert(m->chromosome < perm->n_chrom);
      markers_on_chrs[m->chromosome]++;
    }
  }
  for(i=0; i<perm->n_chrom; i++){
    n_inv += (markers_on_chrs[i]+1)*markers_on_chrs[i];
  }
  n_inv /= 2;
  return n_inv;
}


int
get_n_trans(Permutation* perm)
{
    int NpM = perm->n_mark + perm->n_chrom;
    int NpMchoose2 = NpM*(NpM-1)/2;
    
  return 2*(NpMchoose2 - get_n_inv(perm));
} // end get_n_trans


int
get_n_inv_max(Permutation* perm)
{
        // returns max number of inversions
        // which occurs for case of all markers on 1 chromosome
    int N = perm->n_mark;

    return (N+1)*N/2;
}

int
get_n_inv_min(Permutation* perm)
{
        // returns min number of inversions
        // which occurs for case of markers spread as evenly as possible over chromosomes
        // i.e. N mod M chromosomes with N/M markers, the rest with N/M + 1 markers
    int N = perm->n_mark, M = perm->n_chrom;
    int mpc = N/M;
    int NmodM = N % M;   

    return (M - NmodM)*(mpc+1)*mpc/2 + NmodM*(mpc+2)*(mpc+1)/2;
}

int
get_n_trans_max(Permutation* perm)
{
    int N = perm->n_mark, M = perm->n_chrom;
    int NplusMchoose2 = (N+M)*(N+M-1)/2;

    return (NplusMchoose2 - get_n_inv_min(perm))*2;
}

int
get_n_trans_min(Permutation* perm)
{
    int N = perm->n_mark, M = perm->n_chrom;
    int NplusMchoose2 = (N+M)*(N+M-1)/2;

    return (NplusMchoose2 - get_n_inv_max(perm))*2;
}



// checking functions

int
are_perms_equivalent_both(Cycle_decomposition* the_cd, Permutation* perm1, Permutation* perm2)
{
    int result1 = are_perms_equivalent(perm1, perm2);
    int result2 = are_perms_equivalent_fast(the_cd, perm1, perm2);
    if(result1 != result2) {
        printf("in are_perms_equivalent. old, new ways don't agree: %i %i \n", result1, result2);
        getchar();
    }
    return result1;
        // return are_perms_equivalent_alt(perm1, perm2);
}

int
are_perms_equivalent_fast(Cycle_decomposition* the_cd, Permutation* perm1, Permutation* perm2)
{
        // fastest way if already have cycle decomposition
    int c_int_targ = perm2->n_mark - perm2->n_chrom_nonempty;
   
    if(the_cd->n_int_cycles != c_int_targ) return FALSE;
    if(perm1->n_chrom_nonempty != perm2->n_chrom_nonempty) return FALSE;
  
    return TRUE;  
} // end are_perms_equivalent1

int are_perms_equivalent(const Permutation* perm1, const Permutation* perm2)
{
        // faster way to test - can be further improved by factor of 2 I believe
    int N = perm1->n_mark;
    int i;
    Marker_end* me1, * me2;

    if(perm1 == NULL || perm2 == NULL) return FALSE; 
    if(perm1->n_chrom_nonempty != perm2->n_chrom_nonempty) return FALSE;
    for(i=1, me1 = perm1->p+1, me2 = perm2->p+1; i<=2*N; i++, me1++, me2++){ // loop over all the marker_ends which are not on the ends of chromosomes
        assert(!me1->is_end_marker);
        assert(!me2->is_end_marker);
        if(!me1->black->is_end_marker){
            if(me2->black->is_end_marker || (me2->black->markend_numb != me1->black->markend_numb)) return FALSE;
        }
        else{
            if(!me2->black->is_end_marker) return FALSE;
        }
    }
    return TRUE;
} // end of are_perms_equivalent

int compare_interior_edges(Permutation* perm1, Permutation* perm2, int* N1, int* N2)
{
        // count the number of edges joining two internal (non_end) markers for each perm,
        // and returns the number of these edges common to both perms
    int N = perm1->n_mark;
    int i;
    Marker_end* me1, * me2;
    int n1 = 0, n2 = 0, nboth = 0;
    
    for(i=1, me1 = perm1->p+1, me2 = perm2->p+1; i<=2*N; i++, me1++, me2++){ // loop over all the marker_ends which are not on the ends of chromosomes
        assert(!me1->is_end_marker);
        assert(!me2->is_end_marker);
        if(!me1->black->is_end_marker){ // neighboring int marker_ends
            n1++;
            if(!me2->black->is_end_marker){
                n2++;
                if(me2->black->markend_numb == me1->black->markend_numb){
                    nboth++;
                }
            }    
        }
        else{
            if(!me2->black->is_end_marker) n2++;
        }
    }
    assert(n1%2 == 0); assert(n2%2 == 0); assert(nboth%2 == 0);
    *N1 = n1/2; *N2 = n2/2;
    return nboth/2;
} // end of are_perms_equivalent_alt

/* int are_perms_equivalent_alt(Permutation* perm1, Permutation* perm2) */
/* { */
/*     int N = perm1->n_mark; */
/*     int i; */
/*     Marker_end* me1, * me2; */
    

/*     for(i=1, me1 = perm1->p+1, me2 = perm2->p+1; i<=2*N; i++, me1++, me2++){ // loop over all the marker_ends which are not on the ends of chromosomes */
/*         assert(!me1->is_end_marker); */
/*         assert(!me2->is_end_marker); */
/*         if(!me1->black->is_end_marker){ */
/*             if(me2->black->is_end_marker || (me2->black->markend_numb != me1->black->markend_numb)) return FALSE; */
/*         } */
/*         else{ */
/*             if(!me2->black->is_end_marker) return FALSE; */
/*         } */
/*     } */
/*     return TRUE; */
/* } // end of are_perms_equivalent_alt */

            

int  are_perms_equal(const Permutation* p1, const Permutation* p2)
{
  unsigned int result = 0, factor = 1, all_ok = 0;
  int i, n_checks = 6;
  Marker_end* the_end1; Marker_end* the_end2;
  char same[] = {1,1,1,1, 1,1,1,1};

  if (p1->n_mark != p2->n_mark) same[0] = FALSE;

  for (i=0; i<= 2*p1->n_mark + 1; i++){
    the_end1 = p1->p+i; the_end2 = p2->p+i;
    if (the_end1->markend_numb != the_end2->markend_numb) same[1] = FALSE;
    if (the_end1->position != the_end2->position) same[2] = FALSE;
    if ((the_end1->black - the_end1) != (the_end2->black - the_end2)) same[3] = FALSE;
    if (i>=1 && i <= 2*p1->n_mark) {
      if ((the_end1->other - the_end1) != (the_end2->other - the_end2))  same[4] = FALSE;
    }
    else {
      if (the_end1->other != NULL || the_end2->other != NULL)  same[5] = FALSE;
    }
  }
    
    for (i=0; i<n_checks; i++, factor*=2){
      result += factor*same[i];
      all_ok += factor;
    }
        //if (result != all_ok) {fprintf(stdout, "are_perms_equal result: %o. %o -> perms equal.\n", result, all_ok);  process_error();  }
    return (result == all_ok);
} // end of function are_perms_equal

int check_markend_numbs(const Permutation* perm)
{
	int i;
   int result = TRUE; 
  
	for (i = 0; i < 2*(perm->n_mark + perm->n_chrom); i++){
		if (perm->p[i].markend_numb != i) {result = FALSE;}
	}
	return result;
} // end of function check_markend_numbs

int check_markend_is_end(const Permutation* perm)
{ // check the is_end_marker fields

   Marker_end* end = perm->p, * black;
   
	do{
       black = end->black;
       if(end->chromosome != black->chromosome){
               // look at neighboring pairs. both end markers <-> different chromosome numbers
           if (!(end->is_end_marker && black->is_end_marker)){return FALSE;}
       }
       
	}while((end = black->other) != NULL);
	return TRUE;
} // end of function check_markend_is_end

int check_n_chrom_nonempty(const Permutation* perm)
{ // check the is_end_marker fields
	int n_empty = 0;
   Marker_end* end = perm->p, * black;
   
   
	do{
       black = end->black;
       if(end->is_end_marker && black->is_end_marker){
           n_empty++;
       }
       
	}while((end = black->other) != NULL);
	return (perm->n_chrom == n_empty + perm->n_chrom_nonempty);
} // end of function check_n_chrom_nonempty

int check_markend_positions(const Permutation* perm)
{
  Marker_end* p = perm->p;
  int count = 0;

 while (p != NULL ){
         //  printf("in check_markend_positions. count, p->position, p->markend_number: %i %i  %i \n", count, p->position, p->markend_numb);
    if (p->position != count) return FALSE;
    count++;
    p = p->black;
        //  printf("in check_markend_positions. count, p->position, p->markend_number: %i %i  %i\n", count, p->position, p->markend_numb);
    if (p->position != count) return FALSE;
    count++;
    p = p->other;
  }
  return TRUE;
} // end of function check_markend_positions

int check_markend_black_ptrs(const Permutation* perm)
{
   Marker_end* end;
   int i, pos, black_pos;

   for (i=0, end = perm->p; i<=2*(perm->n_mark + perm->n_chrom)-1; i++, end++){
      pos = end->position;
      black_pos = end->black->position;
      if (pos + 1 - 2*(pos % 2) != black_pos) return FALSE;
      if (end->black->black != end) return FALSE;
   }
   return TRUE;
} // end of function check_markend_black_ptrs

int check_markend_dists(const Permutation* perm)
{
    Marker_end* end, * bend;
    int i, pos, black_pos;
    int retval = TRUE;
        //  int logdiff, logdiffhist[10] = {0};

    for (i=0, end = perm->p; i<=2*(perm->n_mark + perm->n_chrom)-1; i++, end++){
        bend = end->black;
        if(!real_close2(end->d_black, bend->d_black, 5.0e-12, 1.0)){
            printf("end->d_black, end->black->d_black: %g %g \n", end->d_black, end->black->d_black); getchar();
            retval = FALSE;
        }
      /*   logdiff =  (int)(log10(fabs(end->d_black - fabs(end->distance - bend->distance))) + 15.0); */
/*         logdiff = logdiff < 0? 0: logdiff; */
/*         logdiff = logdiff > 9? 9: logdiff; */
/*         logdiffhist[logdiff]++; */
        
        if(!real_close2(end->d_black, fabs(end->distance - bend->distance), 5.0e-12, 1.0)){
            printf("end->d_black, end->distance, end->black->distance: %g %g %g \n",end->d_black, end->distance, end->black->distance); getchar();
            retval = FALSE;
        }
    }
        //  for(i=0; i<10; i++){ printf(" %4i", logdiffhist[i]); }printf("\n");
    
    if(!retval)  print_perm(stdout, perm);
    return retval;
} // end of function check_markend_dists

int check_markend_dists2(const Permutation* perm)
{
    Marker_end* end, * bend;
    int i, pos, black_pos;
    int retval = TRUE;
    double dist = 0.0;

    end = perm->p;
    while(end != NULL){
        bend = end->black;
        if(!real_close2(end->d_black, bend->d_black, 5.0e-12, 1.0)){
            printf("end->d_black, end->black->d_black: %g %g \n", end->d_black, end->black->d_black);    
            retval = FALSE;
        }
        dist += end->d_black;
        if(!real_close2(bend->distance, dist, 1.0e-11, 1.0)){
            printf("bend->distance, dist: %g %g \n", bend->distance, dist);
            retval = FALSE;
        }
        end = bend->other;
    }
   
    if(!retval)  print_perm(stdout, perm);
    return retval;
} // end of function check_markend_dists2


int check_markend_other_ptrs(const Permutation* perm)
{
    Marker_end* end;
    int i, pos, other_pos;
    int result = TRUE;

    for (i=0, end = perm->p; i<=2*(perm->n_mark + perm->n_chrom)-1; i++, end++){
        pos = end->position;
        if(pos == 0 || pos == 2*(perm->n_mark + perm->n_chrom)-1){ // leftmost or rightmost
            result = (result && (end->other == NULL));
        }
        else{
            other_pos = end->other->position;
            if (pos + 1 - 2*((pos+1) % 2) != other_pos) result = FALSE;
            if (end->other->other != end) result = FALSE;
        }
    }
        //if (perm->p->other != NULL) return FALSE;
        // if (perm->p[2*(perm->n_mark)+1].other != NULL) return FALSE;
    return result;
} // end of function check_markend_other_ptrs

int check_perm_selfconsistent(const Permutation* perm)
{
    int result = TRUE;
    result = result && check_markend_numbs(perm); if(!result) printf("after check markend numbs: result: %i \n", result);
    result = result && check_markend_is_end(perm);  if(!result)  printf("after check markend is end: result: %i \n", result);
    result = result && check_n_chrom_nonempty(perm);  if(!result)  printf("after check n_chrom_nonempty: result: %i \n", result);
    result = result && check_markend_positions(perm); if(!result)  printf("after check markend positions: result: %i \n", result); 
    result = result && check_markend_black_ptrs(perm); if(!result)  printf("after check markend black ptrs: result: %i \n", result); 
    result = result && check_markend_other_ptrs(perm);  if(!result)  printf("after check markend other pts: result: %i \n", result);
    result = result && check_markend_dists2(perm); if (!result)  printf("after check markend dists2: result: %i \n", result);
    result = result && check_markend_dists(perm); if (!result)  printf("after check markend dists: result: %i \n", result);

    return result;

} // end of function check_perm_selfconsistent

int check_directions(const Permutation* perm, Cycle_decomposition* the_cd)
{
        //  Cycle* the_cycle = the_cd->cycles;
    Cycle_element* ce;
    int i, j=-100, d1, d2, pos1, pos2, result = TRUE;
    int NpM = perm->n_mark + perm->n_chrom;


   /*  for(i=0; i<the_cd->n_cycles; i++, the_cycle++){ */
/*         for (j=0, ce=&(the_cd->elements[the_cycle->start_elem_index]); j<the_cycle->n_black; j++, ce++){ */
    for(i=0, ce = the_cd->elements; i < NpM; i++, ce++){
            d1 = ce->direction;
            pos1 = perm->p[ce->n1].position;
            pos2 = perm->p[ce->n2].position;
            d2 = (pos1 > pos2);
            if ((d1 != d2)|| (pos1/2 != pos2/2)){
                fprintf(stdout, "In check_directions.n_cycles: %i  cycle #: %i; black edge #: %i, n1, n2: %i %i dir1, dir2: %i %i , pos1, pos2 %i %i \n",
                       the_cd->n_cycles, i, j,ce->n1, ce->n2, d1, d2, pos1, pos2);
               
                print_perm(stdout, perm);
                print_mends(perm);
                result = FALSE;
                process_error();
            }
    }
    return result;
} // end of function check_directions


// allocation, initialization

Permutation* multipermalloc(int n_chromosomes, int n_markers)  // allocates memory for Permutation
{ // including for the array of marker ends - multiple chromosome version
    Permutation* perm;

    perm = (Permutation*)chcalloc(1, sizeof(Permutation)); n_perm_alloc++;
    perm->p = (Marker_end*)chcalloc(2*(n_markers+n_chromosomes), sizeof(Marker_end)); n_p_alloc++;
    perm->n_mark = n_markers; // number of real markers
    perm->n_chrom = n_chromosomes;
    return perm;
} // end of function multipermalloc

void permfree(Permutation** perm)
{
    free((*perm)->p);n_p_alloc--;
    free(*perm); n_perm_alloc--;
    *perm = NULL;    
}  // end of function permfree

Permutation*
init_genome_to_random(int N_chromosomes, int N_markers, int N_randomize)
{
        // just generate some permutation with specified numbers of chromosomes, markers
        // then can randomize by inversions/translocations if needed

    Permutation* perm;
    int markers_per_chrom = N_markers/N_chromosomes;
    int markers_on_this_chromosome;
    int i;
    FILE* fptmp = fopen("tempout", "w"); //
    int i_chrom = 0;

    
    fprintf(fptmp, "Total_N_markers:  %i\n", N_markers);
    fprintf(fptmp, "N_genomes: 1\n");
    fprintf(fptmp, "Genome: 0\n");
    fprintf(fptmp, "Distances? 0\n");
    fprintf(fptmp, "N_chromosomes: %i \n", N_chromosomes);
    for(i=0; i<N_markers; i++){
        if(i % markers_per_chrom == 0 && i_chrom < N_chromosomes){
            i_chrom++;
            markers_on_this_chromosome = markers_per_chrom;
            if(i_chrom == N_chromosomes)
                markers_on_this_chromosome += N_markers - N_chromosomes*markers_per_chrom; 
            fprintf(fptmp, "Chromosome: %i\n", i_chrom);
            fprintf(fptmp, "N_markers: %i Markers:\n", markers_on_this_chromosome);
        }
        fprintf(fptmp, "%i \n", i);
    }
    fclose(fptmp);
    fptmp = fopen("tempout", "r");
    
    perm = multipermalloc(N_chromosomes, N_markers);
        //   printf("in init_genome_to_random. before get_genome_from_file \n");
    get_genome_from_file(fptmp, &perm, FALSE);
        //   printf("in init_genome_to_random. after get_genome_from_file \n"); getchar();
    fclose(fptmp);

    set_Ai_At(perm); // get chromosome lengths, inv/trans areas

    assert(check_markend_dists(perm));
    assert(check_markend_dists2(perm));

   // randomize the marker order
        //  printf("in init_genome_to_random. before marker order randomizing loop \n");
        //  print_perm_markers(stdout, perm);
    
    for(i=0; i<N_randomize; i++){ simple_rand_reverse(perm, 1.0); }
        //  printf("in init_genome_to_random. after marker order randomizing loop \n");
        // print_perm_markers(stdout, perm);
    return perm;
} // end of init_genome_to_random

Permutation*
init_genome(int N_chromosomes, int* N_markers, double* chromosome_length,
            Marker perm_array[][MAX_N_MARKERS_PER_CHROMOSOME], double L_g)
{
        // like perminit, but for multiple chromosomes
        // L_g is the length of genome, use UNKNOWN if distance info not available
    Permutation* perm;
    int i, j;
    Marker_end* m;
    Marker* M;
    int n_nonempty = 0, n_empty=0;
    int N_markers_total = 0;

    for(i=0; i<N_chromosomes; i++){
        N_markers_total += N_markers[i];
        if(N_markers[i] > 0) {n_nonempty++; /* printf("chromosome %i is non_empty\n", i); */ }
        else if(N_markers[i] == 0) {n_empty++; /* printf("chromosome %i is empty\n", i); */ }
        else printf("in init genome, negative # of markers on chromosome\n");
            //  perm->n_markers_on_chromosome[i] = 0;
    }
        //   printf("in init_genome, N_chromosomes, total_N_markers, # nonempty, empty chromosomes: %i %i %i %i\n",
        //        N_chromosomes, N_markers_total, n_nonempty, n_empty);
    perm = multipermalloc(N_chromosomes, N_markers_total);

    
    perm->marker_names = (char**)malloc(sizeof(char*)*perm->n_mark);
    for(i=0; i<perm->n_mark; i++){
            //  printf("allocating marker_name. i: %i \n", i);
        perm->marker_names[i] = (char*)malloc(sizeof(char)*MAXMARKERNAMELENGTH);
    }
    for(i=0; i<perm->n_chrom; i++){
        for(j=0; j<N_markers[i]; j++){
                //   printf("i: %i  j: %i \n", i, j);
            M = perm_array[i]+j;
                //  printf("namesize, name, number:%i_%s_%i_\n", strlen(M->marker_name), M->marker_name, M->marker_number);
                //  strncpy(perm->marker_names[M->marker_number-1], M->marker_name, MAXMARKERNAMELENGTH-1);
            strcpy(perm->marker_names[M->marker_number-1], M->marker_name);
                //  printf("after strncpy \n");
        }
    }

    
    for(i=0; i<N_chromosomes; i++){
        perm->n_markers_on_chromosome[i] = 0;
    }
    
    init_to_given_multi_perm(N_markers, chromosome_length, perm_array, perm);

        //  printf("perm->n_inv: %i \n", perm->n_inv); // hasn't been defined at this point
    
    m = perm->p; m++;
    for(i=1; i<=perm->n_mark; i++, m += 2){
            //   printf("i, chromosome: %i %i \n", i, m->chromosome);
        perm->n_markers_on_chromosome[m->chromosome]++;
    }
    perm->n_inv = 0;
    for(i=0; i<N_chromosomes; i++){
        long int nmoc = perm->n_markers_on_chromosome[i];
        perm->n_inv += ((nmoc+1)*nmoc)/2;
            //  printf("chromosome: %i  ; markers: %i   nmoc+1choose2: %i n_invsofar:  %i\n", i, (int)nmoc,((nmoc+1)*nmoc)/2, perm->n_inv);
    }
        //  printf("n_inv's %i %i \n", perm->n_inv, get_n_inv(perm));
    perm->n_edge = perm->n_mark + perm->n_chrom;
    set_Ai_At(perm);
        //  printf("perm A_i, A_t: %g %g \n", perm->A_i, perm->A_t); 
    assert(perm->n_inv == get_n_inv(perm));
    perm->n_trans = get_n_trans(perm);
    perm->n_chrom_nonempty = n_nonempty;
    if(L_g == 0.0) L_g = N_markers_total + N_chromosomes; // N_markers_total is number of markers
    perm->L_g = L_g; 
        //  printf("in init_genome: n_chrom, n_chrom_nonempty: %i %i \n", perm->n_chrom, perm->n_chrom_nonempty);
    assert(check_perm_selfconsistent(perm));
    return perm;
} // end init_genome



void init_to_given_multi_perm(int* N_markers, double* chromosome_length, Marker in[][MAX_N_MARKERS_PER_CHROMOSOME], Permutation* out)
{
        // n is the number of signed markers 
        // in is array of Markers representing the genome, in[i][j] is jth Marker on ith chromosome
        // out is pointer to  Permutation structure for same permutation
        // for now, assume that out points to an existing structure (i.e. the memory has been allocated)
    int i, j;
    int p_temp[2*NMARKERMAX+2];
    int pos1, pos2, mend1, mend2, mn;
    int me1, me2, b1, b2;
    int N = out->n_mark, M = out->n_chrom;
    int n, msf = 0, l_pos, r_pos, l_mend, r_mend;
    double chromosome_end_distance[MAX_N_CHROMOSOMES+1];

  

                                       
    chromosome_end_distance[0] = 0.0;
    for(j=1; j<=M; j++){
        chromosome_end_distance[j] = chromosome_end_distance[j-1] + chromosome_length[j-1];
            //  printf("j-1: %i chromosome_length[j-1]: %g \n", j-1, chromosome_length[j-1]);
            //  printf("j: %i chromosome_end_distance[j]: %g \n", j, chromosome_end_distance[j]);
    }
    
    
    for(j=0; j<M; j++){ // loop over chromosomes
        n = N_markers[j]; // number of markers on chromosome j
            //   printf("# of markers on chromosome: %i is %i \n", j, n);
        for (i=0; i<n; i++){
            mn = abs(in[j][i].marker_number);   // marker number (before adding in chromosome end markers)    
            pos1 = (i+msf+j)*2+1;
            pos2 = pos1+1;
            mend2 = 2*mn; //   1,2,3 -> 1,3,5   real marker end numbers run from 1 to 2*n_mark
            mend1 = mend2-1; //      -> 2,4,6

            if (in[j][i].marker_number > 0){
                out->p[mend1].position = pos1;
                out->p[mend2].position = pos2;
                p_temp[pos1] = mend1; // gives marker end number of position pos1
                p_temp[pos2] = mend2;
            }
            else if (in[j][i].marker_number < 0){
                out->p[mend1].position = pos2;
                out->p[mend2].position = pos1;
                
                p_temp[pos1] = mend2;
                p_temp[pos2] = mend1;
            }
            else {
                fprintf(stdout, "In init_to_given_multi_perm; in[j][i] = 0; shouldn't happen. i, j: %i %i \n", i, j);
                process_error();
                return;
            }
            out->p[mend1].distance = out->p[mend2].distance = in[j][i].marker_distance;
            
            out->p[mend1].markend_numb = mend1;
            out->p[mend2].markend_numb = mend2;
                // printf("mend1, mend2: %i %i \n", mend1, mend2);

            out->p[mend1].other = out->p+mend2;
            out->p[mend2].other = out->p+mend1;

            out->p[mend1].is_end_marker = FALSE;
            out->p[mend2].is_end_marker = FALSE;

            out->p[mend1].chromosome = j;
            out->p[mend2].chromosome = j;
            
        } // end of loop over markers in chromosome

        l_pos = 2*(msf+j); //  position of left end of this chromosome
        r_pos = 2*(msf+n+j)+1;  //  position of right end of this chromosome

        l_mend = (j==0)? 0:2*(N+j); // marker end numbers of marker_ends on left ends of chromosomes
        r_mend = 2*(N+j)+1;   // marker end numbers of marker_ends on right ends of chromosomes
        

        p_temp[l_pos] = l_mend;
        p_temp[r_pos] = r_mend;

        out->p[l_mend].markend_numb = l_mend;
        out->p[r_mend].markend_numb = r_mend; //printf("l_mend, r_mend: %i %i \n", l_mend, r_mend);
        out->p[l_mend].position = l_pos;
        out->p[r_mend].position = r_pos;
        out->p[l_mend].distance = chromosome_end_distance[j];
        out->p[r_mend].distance = chromosome_end_distance[j+1];

       /*  out->p[l_mend].chromosome = M; */
/*         out->p[r_mend].chromosome = M; */

         out->p[l_mend].chromosome = j;
        out->p[r_mend].chromosome = j; 
    
        if(j==0){ out->p[l_mend].other = NULL;}
        else{ out->p[l_mend].other = out->p+(l_mend-1);}
        if(j==M-1){ out->p[r_mend].other = NULL;}
        else{ out->p[r_mend].other = out->p+(r_mend+1);} 
   
        out->p[l_mend].is_end_marker = TRUE;
        out->p[r_mend].is_end_marker = TRUE;
    
        // do black pointers
        out->p[l_mend].black = out->p + p_temp[l_pos+1];
        out->p[r_mend].black = out->p + p_temp[r_pos-1];

            // do d_black fields (distances to adjacent markers)
        out->p[l_mend].d_black = fabs(out->p[p_temp[l_pos+1]].distance - out->p[l_mend].distance);
        out->p[r_mend].d_black = fabs(out->p[r_mend].distance - out->p[p_temp[r_pos-1]].distance);
    
    
        for (i=2*(msf+j)+1; i<2*(msf+j+n); i+=2){
            me1 = p_temp[i];
            me2 = p_temp[i+1];
            b1 = p_temp[i-1];
            b2 = p_temp[i+2];
            out->p[me1].black = out->p + b1;
            out->p[me2].black = out->p + b2;
                // get d_black fields for these
            out->p[me1].d_black = fabs(out->p[b1].distance - out->p[me1].distance);
            out->p[me2].d_black = fabs(out->p[b2].distance - out->p[me2].distance);      
        }

        msf += n; // number of markers on chromosomes done so far
    }
       
        //  printf("returning from init_to_given_multi_perm\n");
    
}  // end of function init_to_given_multi_perm


// print info
void
print_markers(const Permutation* perm)
{
    int i;
    Marker_end* e;
// print in position order ...
    for(i=0, e=perm->p; i<(perm->n_mark + perm->n_chrom)-1; i++){
        printf(" %6i %6i %6i %6i \n", e->markend_numb, e->position, e->chromosome, e->is_end_marker);
        e = e->black;
         printf(" %6i %6i %6i %6i \n", e->markend_numb, e->position, e->chromosome, e->is_end_marker);
         e = e->other;
    }
     printf(" %6i %6i %6i %6i \n", e->markend_numb, e->position, e->chromosome, e->is_end_marker);
        e = e->black;
         printf(" %6i %6i %6i %6i \n\n", e->markend_numb, e->position, e->chromosome, e->is_end_marker);
}


void
print_edges(Cycle_decomposition* the_cd, const int n_edges)
{
    int i;
    Cycle_element* e;

    for(i=0, e=the_cd->elements; i<n_edges; i++, e++){
        printf(" %6i %6i %6i %6i %6i %6i\n", e->n1, e->n2, e->position, e->direction, e->chromosome, e->eds);
    }printf("\n");
}

void print_perm_markers(FILE* fp, const Permutation* perm)
{
        // prints markers on chromosomes, using marker names

    Marker_end* m = (perm->p)->black; // now m points LH end of first real marker
    Marker_end* other = m->other;// now m points RH end of first real marker

    fprintf(fp, "| ");
    while(other != NULL){
        if(m->is_end_marker) fprintf(fp, " |");
        if(other->is_end_marker) fprintf(fp, "| ");
        
        else if(m->position > other->position){ fprintf(fp, "(%s) ", perm->marker_names[(m->markend_numb-1)/2]); }
        else{ fprintf(fp, "%s ", perm->marker_names[(m->markend_numb-1)/2]); }
        m = other->black; other = m->other;
    } fprintf(fp, " |\n");
}
    


void print_perm(FILE* fp, const Permutation* perm)
{ // prints the marker end numbers in position order
  Marker_end the_end;
//  int Nrecomb;
      
  the_end = perm->p[0];
  
   fprintf(fp, "                      |%5i   %12g %12g\n", the_end.markend_numb, the_end.distance, the_end.d_black); //fprintf(fp, "            %g  %i\n", the_end.distance, the_end.chromosome);
  the_end = *(the_end.black);
  while(the_end.other != NULL){
    fprintf(fp, "%12g    %5i", the_end.d_black, the_end.markend_numb);
   
    if(the_end.is_end_marker) fprintf(fp, " |"); // to separate chromosomes
    else fprintf(fp, " :"); 
    
    the_end = *(the_end.other);

    fprintf(fp, "%5i   %12g %12g\n", the_end.markend_numb, the_end.distance, the_end.d_black); //fprintf(fp, "   %g  %i\n", the_end.distance, the_end.chromosome);
        // Nrecomb = (the_end.is_end_marker || the_end.black->is_end_marker)? -1: the_end.Nrecomb;
        // fprintf(fp, "%5i   %12i %12g\n", the_end.markend_numb, Nrecomb, the_end.d_black); //fprintf(fp, "   %g  %i\n", the_end.distance, the_end.chromosome);

    the_end = *(the_end.black);
  }
  fprintf(fp, "%12g    %5i |        %12g", the_end.d_black, the_end.markend_numb, the_end.distance);  // fprintf(fp, "            %g  %i\n", the_end.distance, the_end.chromosome);
  fprintf(fp, "\n");
//  getchar();
} // end of function print_perm

void print_perm_brief(FILE* fp, const Permutation* perm)
{ // prints the marker numbers in position order, not marker ends. chroms separated by |
  Marker_end* the_end;
      
  the_end = perm->p;
  
   fprintf(fp, "| "); 
  the_end = the_end->black;
  while(the_end->other != NULL){
    if(the_end->is_end_marker) fprintf(fp, "| "); // to separate chromosomes
    else  fprintf(fp, "%i ", (the_end->markend_numb+1)/2);
    the_end = the_end->other->black;
  }
  fprintf(fp, "|");
  fprintf(fp, "\n");
//  getchar();
} // end of function print_perm_brief

void print_perm2(FILE* fp, const Permutation* perm)
{ // prints the marker end numbers in position order
  Marker_end the_end;
//  int Nrecomb;
      
  the_end = perm->p[0];
  
   fprintf(fp, "                      |%5i   %12g %12g\n", the_end.markend_numb, the_end.distance, the_end.d_black); //fprintf(fp, "            %g  %i\n", the_end.distance, the_end.chromosome);
  the_end = *(the_end.black);
  while(the_end.other != NULL){
    fprintf(fp, "%12g    %12s", the_end.d_black, (the_end.is_end_marker)? "chrom_end": perm->marker_names[(the_end.markend_numb-1)/2]);
   
    if(the_end.is_end_marker) fprintf(fp, " |"); // to separate chromosomes
    else fprintf(fp, " :"); 
    
    the_end = *(the_end.other);

    fprintf(fp, "%12s   %12g %12g\n", (the_end.is_end_marker)?
            "chrom_end": perm->marker_names[(the_end.markend_numb-1)/2], the_end.distance, the_end.d_black); //fprintf(fp, "   %g  %i\n", the_end.distance, the_end.chromosome);
        // Nrecomb = (the_end.is_end_marker || the_end.black->is_end_marker)? -1: the_end.Nrecomb;
        // fprintf(fp, "%5i   %12i %12g\n", the_end.markend_numb, Nrecomb, the_end.d_black); //fprintf(fp, "   %g  %i\n", the_end.distance, the_end.chromosome);

    the_end = *(the_end.black);
  }
  fprintf(fp, "%12g    %12s |        %12g", the_end.d_black, (the_end.is_end_marker)?
          "chrom_end": perm->marker_names[(the_end.markend_numb-1)/2], the_end.distance);  // fprintf(fp, "            %g  %i\n", the_end.distance, the_end.chromosome);
  fprintf(fp, "\n");
//  getchar();
} // end of function print_perm
/* char* marker_end_name(const Permutation* p, Marker_end* m) */
/* { */
/*     return (m->is_end_marker)? "chrom_end": p->marker_names[(m->markend_numb-1)/2];\ */
/* } */

void print_perm1(FILE* fp, const Permutation* perm, int Nindiv)
{ // prints the marker end numbers in position order
    Marker_end the_end;
    double theta_hat;
      
    the_end = perm->p[0];
    theta_hat = (the_end.is_end_marker || the_end.black->is_end_marker)? -1.0: (double)the_end.Nrecomb/(double)Nindiv;
    fprintf(fp, "                      |%5i   %12g %12g %12g\n", the_end.markend_numb, theta_hat, the_end.theta, the_end.d_black); //fprintf(fp, "            %g  %i\n", the_end.distance, the_end.chromosome);
    the_end = *(the_end.black);
    while(the_end.other != NULL){
        fprintf(fp, "%12g    %5i", the_end.d_black, the_end.markend_numb);
   
        if(the_end.is_end_marker) fprintf(fp, " |"); // to separate chromosomes
        else fprintf(fp, " :"); 
    
        the_end = *(the_end.other);

            //  fprintf(fp, "%5i   %12g %12g\n", the_end.markend_numb, the_end.distance, the_end.d_black); //fprintf(fp, "   %g  %i\n", the_end.distance, the_end.chromosome);
        theta_hat = (the_end.is_end_marker || the_end.black->is_end_marker)? -1: (double)the_end.Nrecomb/(double)Nindiv;
        fprintf(fp, "%5i   %12g %12g %12g\n", the_end.markend_numb, theta_hat, the_end.theta, the_end.d_black); //fprintf(fp, "   %g  %i\n", the_end.distance, the_end.chromosome);

        the_end = *(the_end.black);
    }
    theta_hat = (the_end.is_end_marker || the_end.black->is_end_marker)? -1: (double)the_end.Nrecomb/(double)Nindiv;
    fprintf(fp, "%12g    %5i |        %12g %12g", the_end.d_black, the_end.markend_numb, theta_hat, the_end.theta);  // fprintf(fp, "            %g  %i\n", the_end.distance, the_end.chromosome);
    fprintf(fp, "\n");
//  getchar();
} // end of function print_perm1



void check_edss(Cycle_decomposition* the_cd)
{
    int i;
    int n_eds = the_cd->n_edss;
    int n_inc, n_nc;
    Cycle* the_eds;
    Cycle_element* ce;
    int fb, lb, nb, ng;

    for (i=0; i<n_eds; i++){
        the_eds = the_cd->edss + i;
            //   assert(the_eds->is_end_cycle == TRUE);
            //  printf("i, the_eds: %i, %p \n", i, the_eds);
        ce = the_eds->start_elem_ptr;
            //    printf("ce: %p \n", ce);
            //   printf("i: %i    %i ~ %i", i, ce->n1, ce->n2);
        ce = ce->next;
        while(ce != the_eds->start_elem_ptr && ce->eds == the_eds->start_elem_ptr->eds){
                //   printf(" - %i ~ %i", ce->n1, ce->n2);
            ce = ce->next;
        } // printf("\n");
        fb = the_eds->first_black; lb = the_eds->last_black;
        nb = the_eds->n_black, ng = the_eds->n_gray;
               printf("     nb,ng,fb,lb, is_end: %i %i  %i %i  %i  \n",
                      the_eds->n_black, the_eds->n_gray, the_eds->first_black, the_eds->last_black, the_eds->is_end_cycle);
        assert(!the_eds->is_end_cycle || ((fb == TRUE || fb == FALSE) && (lb == TRUE || lb == FALSE)));
        assert(the_eds->is_end_cycle || (fb == TRUE && lb == FALSE)); // int cycle is bg
        assert(abs(nb - ng) <= 1);
        assert(!the_eds->is_end_cycle || (!(nb == ng) || (fb != lb)));
        assert(!the_eds->is_end_cycle || (!(nb == ng+1) || ((fb == TRUE) && (lb == TRUE))));
        assert(!the_eds->is_end_cycle || (!(nb == ng-1) || ((fb == FALSE) && (lb == FALSE))));
        assert(the_eds->is_end_cycle || (nb == ng));
    } // printf("\n");
   /*  printf("n_eds, n_empty_in_targ, c_i: %i %i %i   %i\n", */
/*            the_cd->n_edss, the_cd->n_empty_chrom_in_target, the_cd->n_int_cycles, */
/*            the_cd->n_edss + the_cd->n_empty_chrom_in_target - the_cd->n_int_cycles); */
    assert(the_cd->n_edss + the_cd->n_empty_chrom_in_target - the_cd->n_int_cycles == 2*the_cd->n_chrom);
}// end of function check_edss:


void print_edss(Cycle_decomposition* the_cd)
{
    int i;
    int n_eds = the_cd->n_edss;
    int n_inc, n_nc;
    Cycle* the_eds;
    Cycle_element* ce;
        //  int fb, lb, nb, ng;

    for (i=0; i<n_eds; i++){
        the_eds = the_cd->edss + i;
            //   assert(the_eds->is_end_cycle == TRUE);
            //  printf("i, the_eds: %i, %p \n", i, the_eds);
        ce = the_eds->start_elem_ptr;
            //    printf("ce: %p \n", ce);
        printf("i: %i    %i ~ %i", i, ce->n1, ce->n2);
        ce = ce->next;
        while(ce != the_eds->start_elem_ptr && ce->eds == the_eds->start_elem_ptr->eds){
            printf(" - %i ~ %i", ce->n1, ce->n2);
            ce = ce->next;
        }  printf("\n");
       /*  fb = the_eds->first_black; lb = the_eds->last_black; */
/*         nb = the_eds->n_black, ng = the_eds->n_gray; */
/*         printf("     nb,ng,fb,lb, is_end: %i %i  %i %i  %i  \n", */
/*                the_eds->n_black, the_eds->n_gray, the_eds->first_black, the_eds->last_black, the_eds->is_end_cycle); */
    /*     assert(!the_eds->is_end_cycle || ((fb == TRUE || fb == FALSE) && (lb == TRUE || lb == FALSE))); */
/*         assert(abs(nb - ng) <= 1); */
/*         assert(!the_eds->is_end_cycle || (!(nb == ng) || (fb != lb))); */
/*         assert(!the_eds->is_end_cycle || (!(nb == ng+1) || ((fb == TRUE) && (lb == TRUE)))); */
/*         assert(!the_eds->is_end_cycle || (!(nb == ng-1) || ((fb == FALSE) && (lb == FALSE)))); */
/*         assert(the_eds->is_end_cycle || (nb == ng)); */
    }printf("\n");
   /*  printf("n_eds, n_empty_in_targ, c_i: %i %i %i   %i\n", */
/*            the_cd->n_edss, the_cd->n_empty_chrom_in_target, the_cd->n_int_cycles, */
/*            the_cd->n_edss + the_cd->n_empty_chrom_in_target - the_cd->n_int_cycles); */
        //  assert(the_cd->n_edss + the_cd->n_empty_chrom_in_target - the_cd->n_int_cycles == 2*the_cd->n_chrom);
}// end of function print_edss:



void print_mends(const Permutation* perm)
{
    int i;
    Marker_end* the_mend = perm->p;

    fprintf(stdout, "%i ", the_mend->markend_numb);
    while (the_mend != NULL){
        the_mend = the_mend->black;
        
        fprintf(stdout, "%i ",the_mend->markend_numb);
        the_mend = the_mend->other;
        if (the_mend != NULL){
            fprintf(stdout, "%i ",the_mend->markend_numb);
        }
    }
    fprintf(stdout, "(number is number, position is position \n");
    
    for(i=0; i<=2*(perm->n_mark)+1; i++){
        the_mend = perm->p+i;
        fprintf(stdout, "%i ", the_mend->position);
    }
    fprintf(stdout, " (position is number, number is position\n");
} // end of function print_mends

int
check_rev_ordered(const Permutation* perm, Reversal rev)
{
        // Permutation* perm = step->perm;
        //  Reversal rev = step->rev;
    int result = TRUE;
    int i;

    if(rev.marker_numbers[0] != UNKNOWN){
    
        if( (perm->p[rev.marker_numbers[0]].position >= perm->p[rev.marker_numbers[1]].position) ||
            (perm->p[rev.marker_numbers[1]].position >= perm->p[rev.marker_numbers[2]].position) ||
            (perm->p[rev.marker_numbers[2]].position >= perm->p[rev.marker_numbers[3]].position) ) {
            result = FALSE;
            printf("In check_rev_ordered. marker_numbers: %i %i %i %i \n", rev.marker_numbers[0], rev.marker_numbers[1],
                   rev.marker_numbers[2], rev.marker_numbers[3]);
            printf("In check_rev_ordered. positions: %i %i %i %i \n",
                   perm->p[rev.marker_numbers[0]].position, perm->p[rev.marker_numbers[1]].position,
                   perm->p[rev.marker_numbers[2]].position, perm->p[rev.marker_numbers[3]].position);
            print_rev(stdout, rev);
            getchar();
        }
        else{
         
        }
    }
    return result;
}



// end of perm_etc.c




