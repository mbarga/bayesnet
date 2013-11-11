/* Random number generator ran2 from Numerical Recipes.
 * This file is completely self-contained.  The function ran2
 * and its associated counter randx need not be used directly.
 * It is more convenient to use only the functions srandinter
 * and randinter, defined below.
 *
 *  Random number functions:
 *
 *      int srandinter(int seed)        Initialize the random sequence
 *                                          (seed = 0 ==> use system clock).
 *                                      Return value of seed actually used.
 *
 *      float randinter(float a,
 *                      float b)        Return a random number uniformly
 *                                      distributed between a and b.
 *
 * Cut and paste this file into your program, or add the functions
 * below to your own library.
 *
 */

/*---------------------------- (cut from here...) ---------------------------*/
#include <stdlib.h>

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
#define EPS 1.e-14
#define RNMX (1.0-EPS)

double ran2(int *idum)
{
    int j;
    int k;
    static int idum2 = 123456789;
    static int iy = 0;
    static int iv[NTAB];
    double temp;

    if (*idum <= 0) {           // *idum < 0 ==> initialize
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);

        for (j = NTAB+7; j >= 0; j--) {
            k = (*idum)/IQ1;
            *idum = IA1*(*idum-k*IQ1) - k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum = IA1*(*idum-k*IQ1) - k*IR1;
    if (*idum < 0) *idum += IM1;

    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;

    j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;

    if ((temp = AM*iy) > RNMX)
        return RNMX;
    else
        return temp;
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

static int randx = 0;           /* copy of random seed (internal use only) */

#include <time.h>

int srandinter(int seed)        /* initialize the random number generator */
{
    if (seed == 0) seed = (int) time(NULL);     /* initialize from the system
                                                   clock if seed = 0 */
    randx = -abs(seed);
    return seed;                /* return seed in case we need to repeat */
}

float randinter(float a, float b)       /* return a random number uniformly
                                           distributed between a and b */
{
    if (randx == 0) srandinter(0);
    return a + (b-a)*((float)ran2(&randx));
}

/*------------------------------ (...to here) ------------------------------*/
