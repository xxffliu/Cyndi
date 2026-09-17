#ifndef RANDOM_H
#define RANDOM_H
#include <math.h>
#include <time.h>
#include <cstdlib>

 
/* v4 (2026): these were declared 'static' in the header, which gives EVERY
   translation unit that includes random.h its own private copy of the
   generator state.  It works today only because randomize()/randomperc()/
   advance_random() all live in random.cpp and therefore all touch random.cpp's
   copy -- the moment any other file referenced oldrand it would be seeding a
   different generator than the one it draws from.  Declare them extern here
   and define them once in random.cpp. */
extern double oldrand[55];                      /* Array of 55 random numbers */
extern double rndx1, rndx2;                     /* used with random normal deviate */
extern int rndcalcflag;                         /* used with random normal deviate */

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
