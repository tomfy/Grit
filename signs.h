// signs.h  declaration of functions in signs.c//Dealing with signs for unsigned case

// macros related to signs
#define FLIPONERANDOM 1
#define FLIPONEUSECD 3
#define FLIPSOMEUSECD 6
#define FLIPUSECDMIX  8 // a mix of 3 and 6 i.e. 50% of time 3, 50% 6
#define DELTAD0_FLIP_PROB 0.5  // 
#define DELTADM1_FLIP_PROB 0.7  // Convergence is fairly insensitive to this in range 0.5 to 0.9
#define GAMMA 1.0 // For choose_signs == FLIPONEUSECD, prob of delta c = -1 is GAMMA*Flip_one_sign_epsilon^2
/* #define DDM1 -1 */
/* #define DD0SEDS 0 */
/* #define DD0DEDS 0 */
/* #define DD1 1 */
#define NOTFLIP -1000

void flip_signs_by_position(Permutation* p1, int* flips);
void flip_signs_by_number(Permutation* p1, int* flips);
double get_flipped_signs(int n, double* flip_probs, int* flipped);
void fix_signs(Path* path, Step* last, int* flipped);
double randomize_signs(int n, double pflip, int* flipped);
double get_flips(int n, double* flip_probs, int* flips);
double get_flips_prob(int n, double* flip_probs, int* flips);
void get_flips_by_number(Permutation* perm, int* flips_by_pos, int* flips_by_number);
void flip_one(Permutation* perm, int* flips);
void get_flip_deltacs(Permutation* perm1, Permutation* perm2, int* deltacs);
void get_flip_probs(Permutation* p1, Permutation* p2, double sign_epsilon, double deltac1_flip_prob, double* flip_prob);
double flip_signs(const Run_info_in* r_in, Permutation* perm1old, Permutation* perm1new, Permutation* perm2, int* flips);
int get_signs_short(Permutation* p1, Permutation* p2, const Run_info_in* r_in);

void
get_flip_deltads(Permutation* perm1, Permutation* perm2, int* deltads, int* seds);
double
choose_one_to_flip(Permutation* perm1, Permutation* perm2, double sign_epsilon, int* flips);
double
prob_one_to_flip(Permutation* perm1, Permutation* perm2, double sign_epsilon, int* flips);
// void fix_rev(Reversal* rev, int* flips_by_number);
void
fix_rev(Permutation* perm, Reversal* rev, int* flips_by_number);
int
get_conserved_flips(Permutation* perm1, Permutation* perm2, int* flips);

// end of signs.h
