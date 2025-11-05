//#include "exts.h"
#include "grit.h"
#include "rngs.h"

float ran0(long* idum);
float ran1(long* idum);
float ran2(long* idum);
float ran3(long* idum);
float ran4(long* idum);
double genrand(void);
void sgenrand(unsigned long seed);

unsigned long rng_long; // external variable 

#define RNG 5 // 0,1,2,3,4 Numerical recipes RNGs ran0, ran1, etc; 5 mersenne twister; -1 -> quick_n_dirty, 
#define RANDMAX RAND_MAX

/* note #undef's at end of file */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long *idum)
{
	long k;
	float ans;

	*idum ^= MASK;
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	ans=AM*(*idum);
	*idum ^= MASK;
       // printf("%i %g %g \n", *idum, AM, ans);
	return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
// end of ran0

/* note #undef's at end of file */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
// end of ran1

/* note #undef's at end of file */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
// end of ran2

#include <stdlib.h>
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=labs(MSEED-labs(*idum));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
// end of ran3

float ran4(long *idum)
{
	void psdes(unsigned long *lword, unsigned long *irword);
	unsigned long irword,itemp,lword;
	static long idums = 0;
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static unsigned long jflone = 0x00004080;
	static unsigned long jflmsk = 0xffff007f;
#else
	static unsigned long jflone = 0x3f800000;
	static unsigned long jflmsk = 0x007fffff;
#endif

	if (*idum < 0) {
		idums = -(*idum);
		*idum=1;
	}
	irword=(*idum);
	lword=idums;
	psdes(&lword,&irword);
	itemp=jflone | (jflmsk & irword);
	++(*idum);
	return (*(float *)&itemp)-1.0;
}
// end of ran4

/* note #undef's at end of file */
#define NITER 4

void psdes(unsigned long *lword, unsigned long *irword)
{
	unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
	static unsigned long c1[NITER]={
		0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	static unsigned long c2[NITER]={
		0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

	for (i=0;i<NITER;i++) {
		ia=(iswap=(*irword)) ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		ib=itmpl*itmpl+ ~(itmph*itmph);
		*irword=(*lword) ^ (((ia = (ib >> 16) |
			((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
		*lword=iswap;
	}
}
#undef NITER
// end of psdes (psdes = pseudo data encryption std, needed for ran4)


// mersenne twister rng

/* A C-program for MT19937: Real number version                */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* Any feedback is very welcome. For any question, comments,       */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email       */
/* matumoto@math.keio.ac.jp                                        */

#include<stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void
sgenrand(unsigned long seed)
//  unsigned long seed; 
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

double /* generating reals */
/* unsigned long */ /* for integer generation */
genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

        // return ( (double)y / (unsigned long)0xffffffff ); //reals
        // return y;   //for integer generation

        // return y*(1.0/4294967295.0); // reals:  [0,1]
    
    return y*(1.0/4294967296.0); /* reals: [0,1) */ 
}


double drand(unsigned long* rng_val)
{
    float x;
   
    if(RNG == -1){ return (double)RANDF()/((double)RANDMAX + 1); } // built in rand function
    else if(RNG == 0) { x = ran0((long*)rng_val); } // ran0-ran4 from numerical recipes
    else if(RNG == 1) { x = ran1((long*)rng_val); }
    else if(RNG == 2) { x = ran2((long*)rng_val); }
    else if(RNG == 3) { x = ran3((long*)rng_val); }
    else if(RNG == 4) { x = ran4((long*)rng_val); }
    else if(RNG == 5){ return genrand();} // mersenne twister
    else {
        printf("in drand, unknown rng option selected\n");
        exit(EXIT_FAILURE);
    }
    return (double)x;
} // end of function drand

void seedrand(unsigned long seed)
{
    long lseed = seed; //signed
    if(RNG == -1){ srand(seed); }
    else{ // call ranx with (ptr to) negative value for x = 1,2,3,4
            // call sgenrand for RNG == 5
        lseed = (lseed >= 0)? lseed: -lseed; // absolute value
        rng_long = -lseed; // now rng_long is <= 0, seed is >= 0
        printf("in seedrand, seed, rng_long = %i %i  \n", seed, rng_long);
        if(RNG == 0) {printf("in seedrand, ran0 returned: %f \n", ran0((long*)&rng_long));}
        else if(RNG == 1) {printf("in seedrand, ran1 returned: %f \n", ran1((long*)&rng_long));}
        else if(RNG == 2) {printf("in seedrand, ran2 returned: %f \n", ran2((long*)&rng_long));}
        else if(RNG == 3) {printf("in seedrand, ran3 returned: %f \n", ran3((long*)&rng_long));}
        else if(RNG == 4) {printf("in seedrand, ran4 returned: %f \n", ran4((long*)&rng_long));}
        else if(RNG == 5) {sgenrand(seed);}
        else {
            printf("in drand, unknown rng option selected\n");
            exit(EXIT_FAILURE);
        }
    }
        //  printf("sizeof float, double: %i %i \n", sizeof(float), sizeof(double));
} // end of function seedrand
