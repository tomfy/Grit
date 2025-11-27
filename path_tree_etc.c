// path_tree_etc.c   path_tree and Histogram stuff

#include "grit.h"
#include "structs.h"
#include "exts.h"
#include "path.h"
#include "path_tree_etc.h"
#include "gr_misc.h"



void
path_tree_node_insert_multi(Path_tree_node** tree_node_ptrptr, State* state, int weight)
{     // insert state into tree which keeps track of how many times each path has
        // occurred (for short paths) and how many times each path length has occurred for long paths.
        // calls itself recursively
       
    int compare, count;
    Path* the_path = state->path;
    Path_tree_node* tree_node = *tree_node_ptrptr;
    
    if(tree_node == NULL) // empty (sub)tree, insert new tree_node at top
    {
      tree_node = (Path_tree_node*)chcalloc(1, sizeof(Path_tree_node)); n_path_tree_node_alloc++;
        tree_node->revseq = NULL;
        tree_node->path = copy_path(the_path);
        tree_node->count = weight;
            // rel_pi only meaningful for r=2 and fixed Lambda (i.e. Modes 0,2,4).
        tree_node->pi_rel = rel_pi(state); // exp(-lambda)lambda^L/L! * (N+M choose 2)^-L * 1/2^L_t
        tree_node->pi_rel_old = tree_node->pi_rel*pow(2.0, (double)the_path->length_f); // to account for lumping together of 2^l_f paths in path_tree
            //   tree_node->pi_rel *= (double)the_path->fission_factor;
        tree_node->fission_factor = the_path->fission_factor;
        tree_node->swap_end_factor = the_path->swap_end_factor;
               tree_node->pi_rel *= the_path->swap_end_factor*the_path->fission_factor;
        tree_node->left = tree_node->right = NULL;
        *tree_node_ptrptr = tree_node;
    }
    else if((compare = compare_paths(the_path, tree_node->path, &count)) == 0){ // increment count of appropriate bin
       
        tree_node->count += weight;
        assert(the_path->length == tree_node->path->length);
            //   assert(the_path->length_t == tree_node->path->length_t);
            // if(!real_close(rel_pi(state)*pow(2.0, (double)the_path->length_f), tree_node->pi_rel));
        if(state->Mode == 0 || state->Mode == 2 || state->Mode == 4 ){
                //  printf("length_f's: %i %i \n", the_path->length_f, tree_node->path->length_f);
            //   printf("rel_pi, pi_rel: %g %g \n", rel_pi(state), tree_node->pi_rel);
                //  printf("rel_pi, pi_rel: %g %g \n", rel_pi(state)*state->path->fission_factor*state->path->swap_end_factor, tree_node->pi_rel);
            assert(the_path->length_f == tree_node->path->length_f);
            assert(real_close(rel_pi(state)*state->path->fission_factor*state->path->swap_end_factor, tree_node->pi_rel, DBL_EPSILON*500.0));
            assert(real_close(state->path->fission_factor, tree_node->fission_factor, DBL_EPSILON*500.0));
            assert(real_close(state->path->swap_end_factor, tree_node->swap_end_factor, DBL_EPSILON*500.0));
            assert(real_close(rel_pi(state)*pow(2.0, (double)the_path->length_f), tree_node->pi_rel_old, DBL_EPSILON*RCFACTOR_1));
        }
    }
    else    // if (compare >0) // insert into right sub-tree -
            // left subtrees never used because compare just tells equal or not, doesn't give > or < info
            // seems to be no way to order paths which is invariant under relabelling chromosome end markers, etc..
            // well, could probably do by defining canonical form for perm, but haven't done it ...
    {
        path_tree_node_insert_multi(&(tree_node->right), state, weight);
    }     
} // end of path_tree_node_insert_multi

void
free_path_tree(Path_tree* tree)
{ 
    free_path_subtree(tree->start);    
    if(tree != NULL) free(tree); n_path_tree_alloc--;
} // end of free_path_tree

void
free_path_subtree(Path_tree_node* tree_node)
{ 
        // Path_tree_node* bin;
    if(tree_node == NULL) return;
    if(tree_node->left != NULL){ // non NULL left subtree - free it
        free_path_subtree(tree_node->left); 
    }
    else if (tree_node->right != NULL){ // non NULL right subtree - free it
        free_path_subtree(tree_node->right); 
    }
    else { // node is a leaf - free it
        free(tree_node->revseq); n_revseq_alloc--;
        free(tree_node); n_path_tree_node_alloc--;
    }
} // end of free_path_subtree


Path_tree*
alloc_path_tree(void)
{ 
    Path_tree* tree;
 
    tree = (Path_tree*)chcalloc(1, sizeof(Path_tree)); n_path_tree_alloc++;
    tree->start = NULL;
    return tree; 
} // end of alloc_path_tree


Path_tree_node*
extract_smallest_from_tree(Path_tree* tree)
{
    Path_tree_node* result;
    
    if(tree->start == NULL){ // tree is empty
        return NULL;
    }
    else{
        result = extract_smallest_from_subtree(&(tree->start));
        return result;        
    }
}  // end of extract_smallest_from_tree 
    
Path_tree_node*
extract_smallest_from_subtree(Path_tree_node** tree_node_ptrptr)
{
    Path_tree_node* tree_node_ptr = *tree_node_ptrptr;
    
    if(tree_node_ptr->left != NULL){  // extract from left subtree
        return extract_smallest_from_subtree(&(tree_node_ptr->left));
    }
    else{ // take this node, promote right subtree (which may be NULL)
        *tree_node_ptrptr = tree_node_ptr->right;
        return tree_node_ptr;
    }
} // end of extract_smallest_from_subtree

void
print_path_tree(FILE* fp, Path_tree* tree, int print_it)
{
        // prints out the contents of the path tree destructively,
        // i.e. the node are removed from the tree which ends up empty
        // if print_it is FALSE then just empty tree, don't print
    int i=1, printed_something = 0;
    Path_tree_node* tree_node;

    while((tree_node = extract_smallest_from_tree(tree)) != NULL){
        if(print_it) {
            printed_something = 1;
            fprintf(fp, "i: %4i  l: %4i  l_t: %4i  l_f: %4i  fission_factor: %4g   counts: %8i pirel, counts/pi_rel: %g %g \n",
                             i, tree_node->path->length,  tree_node->path->length_t, tree_node->path->length_f,
                    tree_node->path->fission_factor, tree_node->count,
                    tree_node->pi_rel, (double)tree_node->count/tree_node->pi_rel);
                //  print_path(fp, tree_node->path); fprintf(fp, "\n");
        }
            // if(print_it) fprintf(fp, "i: %i n_inv: %i counts: %i \n", i, tree_node->revseq[0], tree_node->count);
        i++;
            // free(tree_node->revseq); n_revseq_alloc--;
        free(tree_node); n_path_tree_node_alloc--;
    }
        //if(printed_something)
    fprintf(fp, "\n");
}  // end of print_path_tree

/* Histogram*  */
/* create_histogram(int nbins, double minval, double binwidth) */
/* { */
/*   Histogram* hist; */

/*   hist = (Histogram*)malloc(sizeof(Histogram)); */
/*   hist->nbins = nbins; */
/*   hist->minval = minval; // lower edge of leftmost bin, i.e. smallest value that goes in bin hist->histogram[0] */
/*   hist->binwidth = binwidth; */
/*   hist->sum = 0.0; */
/*   hist->peak = 0.0; */
/*   hist->peak_bin = -1; */
/*   hist->mode = UNKNOWN; */
/*   hist->histogram = (double*)calloc((size_t)nbins, sizeof(double)); */
/*   hist->stats = (Stats*)malloc(sizeof(Stats)); */
/*   initialize_stats(hist->stats); */
/*   return hist; */
/* } */

/* void */
/* cred_int(Histogram* hist, double CL, double* xlo, double* xhi, double* ylo, double* yhi) */
/* { */
/*     int i=0,j=hist->nbins-1; */
/*     double* a=hist->histogram;  */
/*     double prob=0.0, newprob, pi, pj; */
/*     double sum = hist->sum, binwidth = hist->binwidth, tailprob=1.0-CL; */

/*     *xlo = hist->minval; *xhi = *xlo + binwidth*hist->nbins; */
/*     *ylo = *yhi = 0.0; */
    
/*     do{ */
/*         pi = a[i]/sum; */
/*         pj = a[j]/sum; */
/*         if(pi < 0.0 || pj < 0.0) { */
/*             printf("in cred_int, negative probability\n"); */
/*         } */
/*         if (pi < pj){ */
/*             newprob = prob + pi; */
/*             *ylo = pi/binwidth; */
/*             if (newprob > tailprob){ // *xlo will lie in ith bin */
/*                 *xlo = *xlo + (tailprob - prob)/pi*binwidth; printf("xlo, xhi: %g %g \n", *xlo, *xhi); */
/*                     // *ylo = pi/binwidth; *yhi = pj/binwidth; */
/*                 return; */
/*             } */
            
/*             prob = newprob; */
/*             i++; *xlo += binwidth; */
/*         } */
/*         else if (pj <= pi){ */
/*             newprob = prob + pj; */
/*             *yhi = pj/binwidth; */
/*             if (newprob > tailprob){ // *xhi will lie in jth bin */
/*                 *xhi = *xhi - (tailprob - prob)/pj*binwidth; printf("xlo, xhi: %g %g \n", *xlo, *xhi); */
/*                     // *ylo = pi/binwidth; *yhi = pj/binwidth; */
/*                 return; */
/*             } */
/*             prob = newprob; */
/*             j--; *xhi -= binwidth; // printf("xhi: %g \n", *xhi); */
/*         } */
/*     }while(TRUE); */
/* }  // end of cred_int */

/* void */
/* cred_int1(double* x, double* a, int N, double CL, double* xlo, double* xhi, double* ylo, double* yhi) */
/* { */
/*         // gets credible interval with confidence level CL, given arrays x and a */
/*     int i,j; */
/*     double prob=0.0, newprob, pi, pj; */
/*     double sum = 0.0, binwidth = x[1] - x[0], tailprob=1.0-CL; */

/*     *xlo = x[0] - 0.5*binwidth; *xhi = x[N-1] + 0.5*binwidth; */
/*     *ylo = *yhi = 0.0; */

/*     for(i=0; i<N; i++){ */
/*         sum += a[i]; */
/*     } */
        
/*     i=0; j=N-1; */
/*     do{ */
/*         pi = a[i]/sum; */
/*         pj = a[j]/sum; */
/*         if(pi < 0.0 || pj < 0.0) { */
/*             printf("in cred_int, negative probability\n"); */
/*         } */
/*         if (pi < pj){ */
/*             newprob = prob + pi; */
/*             *ylo = pi/binwidth; */
/*             if (newprob > tailprob){ // *xlo will lie in ith bin */
/*                 *xlo = *xlo + (tailprob - prob)/pi*binwidth; printf("xlo, xhi: %g %g \n", *xlo, *xhi); */
/*                     // *ylo = pi/binwidth; *yhi = pj/binwidth; */
/*                 return; */
/*             } */
            
/*             prob = newprob; */
/*             i++; *xlo += binwidth; */
/*         } */
/*         else if (pj <= pi){ */
/*             newprob = prob + pj; */
/*             *yhi = pj/binwidth; */
/*             if (newprob > tailprob){ // *xhi will lie in jth bin */
/*                 *xhi = *xhi - (tailprob - prob)/pj*binwidth; printf("xlo, xhi: %g %g \n", *xlo, *xhi); */
/*                     // *ylo = pi/binwidth; *yhi = pj/binwidth; */
/*                 return; */
/*             } */
/*             prob = newprob; */
/*             j--; *xhi -= binwidth; // printf("xhi: %g \n", *xhi); */
/*         } */
/*     }while(TRUE); */
/* }  // end of cred_int1 */
        
    

/* void */
/* free_histogram(Histogram* hist) */
/* { */
/*     free(hist->stats); */
/*   free(hist->histogram); */
/*   free(hist); */
/* } */

/* void  */
/* insert_in_histogram(Histogram* hist, double value, double weight) */
/* { */
/*   int bin = (int)((value - hist->minval)/(hist->binwidth)); */
/*   if (bin < 0) bin = 0; */
/*   if (bin >= hist->nbins) bin = (hist->nbins)-1; */
/*   hist->histogram[bin] += weight; */
/*   hist->sum += weight; */
/*   if(hist->histogram[bin] > hist->peak){ */
/*       hist->peak = hist->histogram[bin]; */
/*       hist->peak_bin = bin; */
/*       hist->mode = hist->minval + (0.5 + bin)*hist->binwidth; */
/*   } */
/*   insert_in_stats(hist->stats, value, weight); */
/* } */

/* void */
/* print_histogram(FILE* fp, Histogram* hist, int normalize, double A) */
/* { // print a histogram */
/*         // normalize = 0    don't normalize */
/*         // normalize = 1  normalize so sum of bins is A */
/*         // normalize = 2   normalize so peak is A */
/*         // normalize = 3  normalize so sum of bins is A/binwidth (i.e. integrate) */
/*     int i; */
/*     double denom = 1.0; */
/*     double v, max = UNKNOWN; */
/*     double* h; */
/*     int lo_print_bin = 0, hi_print_bin = hist->nbins - 2; */


/*     if(hist->sum > 0.0){ */
/*         for(i=1; i<hist->nbins-1; i++){ */
/*             if(hist->histogram[i] > 0){ lo_print_bin = i-1; break; } */
/*         } */
/*         for(i=hist->nbins-3; i>0; i--){ */
/*             if(hist->histogram[i] > 0){ hi_print_bin = i+1; break; } */
/*         } */
        
/*         if(normalize == 1){ */
/*             if(hist->sum > 0.0)denom = hist->sum/A; */
/*             else { */
/*                 printf("In print_histogram. Can't normalize histogram, sum is <= 0\n"); */
/*                 denom = 1.0; */
/*             } */
/*         } */
/*         else if(normalize == 2){ */
/*             for(i=0, h = hist->histogram; i<hist->nbins; i++, h++){ */
/*                 v = *h; */
/*                 if (v > max) max = v; */
/*             } */
/*             if(max > 0.0){ */
/*                 denom = max/A; */
/*             } */
/*             else  { */
/*                 printf("In print_histogram. Can't normalize histogram, max is <= 0\n"); */
/*                 denom = 1.0; */
/*             } */
/*         } */
/*         else if(normalize == 3){ */
/*             if(hist->sum > 0.0)denom = hist->sum*hist->binwidth/A; */
/*             else { */
/*                 printf("In print_histogram. Can't normalize histogram, sum is <= 0\n"); */
/*                 denom = 1.0; */
/*             } */
/*         } */
      
/*     } */
/*     else{ */
/*         lo_print_bin = 0; hi_print_bin = 1; denom = 1.0; */
/*     } */

  
 
  
/*     for(i=lo_print_bin; i<=hi_print_bin; i++){ // last bin is overflow - don't print */
/*         fprintf(fp, "%g %g %g\n", hist->minval + hist->binwidth*(i + 0.5), hist->histogram[i], hist->histogram[i]/denom); */
/*     } */
/*     fprintf(fp, "\n"); */
/* } */

/* void */
/* print_histograms(FILE* fp, int Nhists, Histogram** the_hists, int normalize, double A) */
/* { // print a histogram */
/*         // normalize = 0    don't normalize */
/*         // normalize = 1  normalize so sum of bins is A */
/*         // normalize = 2   normalize so peak is A */
/*         // normalize = 3 normalize so integral is A */
/*     int i, j; */
/*     double v, max = UNKNOWN; */
/*     double* h; */
/*     Stats the_stats; */
/*     Histogram* hist; */
/*     double denom[NCHAINMAX]; */
/*     int lo_print_bin[NCHAINMAX], hi_print_bin[NCHAINMAX]; */
/*     int lpb, hpb; */
/*     int non_empty = FALSE; //  */

/*     for(j=0; j<Nhists; j++){ */
/*         if(the_hists[j]->sum > 0.0){ non_empty = TRUE; } */
/*     } */
  
/*     if(non_empty){ */
/*         for(j=0; j<Nhists; j++){  */
/*             for(i=1; i<the_hists[j]->nbins-1; i++){ */
/*                 if(the_hists[j]->histogram[i] > 0){ lo_print_bin[j] = i-1; break; } */
/*             } */
/*             for(i=the_hists[j]->nbins-3; i>0; i--){ */
/*                 if(the_hists[j]->histogram[i] > 0){ hi_print_bin[j] = i+1; break; } */
/*             } */
/*         } */
/*         lpb = lo_print_bin[0]; hpb = hi_print_bin[0]; */
/*         for(j=1; j<Nhists; j++){ */
/*             if(lo_print_bin[j] < lpb) lpb = lo_print_bin[j]; */
/*             if(hi_print_bin[j] > hpb) hpb = hi_print_bin[j]; */
/*         } */

      
/*             // normalize all the histograms */
/*         for (j=0; j<Nhists; j++){ */
/*             hist = the_hists[j]; */
/*             if(normalize == 1){ */
/*                 if(hist->sum > 0.0)denom[j] = hist->sum/A; */
/*                 else { */
/*                     printf("In print_histogram. Can't normalize histogram, sum is <= 0\n"); */
/*                     denom[j] = 1.0; */
/*                 } */
/*             } */
/*             else if(normalize == 2){ */
/*                 for(i=0, h = hist->histogram; i<hist->nbins; i++, h++){ */
/*                     v = *h; */
/*                     if (v > max) max = v; */
/*                 } */
/*                 if(max > 0.0){ */
/*                     denom[j] = max/A; */
/*                 } */
/*                 else  { */
/*                     printf("In print_histogram. Can't normalize histogram, max is <= 0\n"); */
/*                     denom[j] = 1.0; */
/*                 } */
/*             } */
/*             else if(normalize == 3){ */
/*                 if(hist->sum > 0.0)denom[j] = hist->sum*hist->binwidth/A; */
/*                 else { */
/*                     printf("In print_histogram. Can't normalize histogram, sum is <= 0\n"); */
/*                     denom[j] = 1.0; */
/*                 } */
/*             } */
/*         } */
/*     } */
/*     else{ */
/*         printf("In print_histograms. All histograms are empty. \n"); */
/*         lpb = 0; hpb = 1; */
/*         for(j=0; j<Nhists; j++){ denom[j] = 1.0; } */
/*     } */

       
/*     for(i=lpb; i<=hpb; i++){ */
/*         initialize_stats(&the_stats); */
/*         for(j=0; j<Nhists; j++){           */
/*             insert_in_stats(&the_stats, the_hists[j]->histogram[i]/denom[j], 1.0); */
/*         } */
/*         fprintf(fp, "%g %g %g ", the_hists[0]->minval +  the_hists[0]->binwidth*(i + 0.5), */
/*                 the_stats.mean, the_stats.stddev_mean); */
     
/*         for(j=0; j<Nhists; j++){ */
/*             fprintf(fp, "%g ", the_hists[j]->histogram[i]/denom[j]); */
/*         }fprintf(fp, "\n"); */
/*     }fprintf(fp, "\n"); */
/* //  fprintf(fp, "%g %g \n\n", the_hists[0]->minval + the_hists[0]->binwidth*(20.0*hist->nbins + 0.5), 0.0); */
/* } // end print_histograms */

/* void */
/* initialize_stats(Stats* stat) */
/* { */
/*     stat->N = 0; stat->sum = 0.0; stat->sumsqr = 0.0; */
/* } */

/* void */
/* insert_in_stats(Stats* stat, double value, double weight) */
/* { */
/*     double mean; */
/*     stat->N += weight; */
/*     stat->sum += weight*value; */
/*     stat->sumsqr += weight*value*value; */
/*     mean = stat->sum/stat->N; */
/*     stat->mean = mean; */
/*     if(stat->N>1){ */
/*     stat->stddev_mean = sqrt((stat->sumsqr/stat->N - mean*mean)/(stat->N-1)); */
/*     stat->stddev = stat->stddev_mean*sqrt(stat->N); */
/*     } */
/*     else{ */
/*         stat->stddev = stat->stddev_mean = UNKNOWN; */
/*     } */
/* } */

/* void print_stats(FILE* fp, Stats* stats) */
/* { */
/*     fprintf(fp, "N: %g; mean: %g; sigma: %g; sigma_mean: %g \n", stats->N, stats->mean, stats->stddev, stats->stddev_mean); */
/* } */

