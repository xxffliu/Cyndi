#ifndef RANDOM_H
#define RANDOM_H
#include <math.h>
#include <time.h>
#include <cstdlib>

 
/* variables are declared static so that they cannot conflict 
with names of   */ 
/* other global variables in other files.  See K&R, p 80, for 
scope of static */
static double oldrand[55];                      /* Array of 55
random numbers */
//static int jrand;                                    /*
//current random number */
static double rndx1, rndx2;                       /* used with random
normal deviate */
static int rndcalcflag;                    /* used with random
normal deviate */

void advance_random(void);


int flip(float prob);

void randomize(float randomseed);
 
double randomnormaldeviate();


double randomperc();


int rnd(int low, int high);


inline double rndreal(float lo ,float hi)
/* real random number between specified limits */
//float lo, hi;
{
    return((randomperc() * (hi - lo)) + lo);
}


void warmup_random(float random_seed);

double noise(double mu ,double sigma);

void initrandomnormaldeviate();

float TimeRandomSeed();
#endif
/*-------------------------------------------------------------*/
