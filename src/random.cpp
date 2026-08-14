#include "../include/random.h"

static int jrand(0);                                    /*
current random number */

void initrandomnormaldeviate()
/* initialization routine for randomnormaldeviate */
{
    rndcalcflag = 1;
}

double randomperc()
/* Fetch a single random number between 0.0 and 1.0 -  */
/* Subtractive Method . See Knuth, D. (1969), v. 2 for */
/* details.Name changed from random() to avoid library */
/* conflicts on some machines                          */
{
    jrand++;
    if(jrand >= 55)
    {
        jrand = 1;
        advance_random();
    }
    return(oldrand[jrand]);
}

 
void advance_random()
/* Create next batch of 55 random numbers */
{
    int j1;
    double new_random;

    for(j1 = 0; j1 < 24; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1+31];
        if(new_random < 0.0) new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    }
    for(j1 = 24; j1 < 55; j1++)
    {
        new_random = oldrand [j1] - oldrand [j1-24];
        if(new_random < 0.0) new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    }
}


int flip(float prob)
/* Flip a biased coin - true if heads */
{
    //float randomperc();

    if(randomperc() <= prob)
        return(1);
    else
        return(0);
}


void randomize(float seed)
/* Get seed number for random and start it up */
{
    int j1;

    for(j1=0; j1<=54; j1++)
		oldrand[j1] = 0.0;
    //jrand=0;

    warmup_random(seed);
    return;
}

double randomnormaldeviate()
/* random normal deviate after ACM algorithm 267 / Box-Muller Method */
{
    //double randomperc();
    double t;

    if(rndcalcflag)
    {
        rndx1 = sqrt(- 2.0*log((double) randomperc()));
        t = 6.2831853072 * (double) randomperc();
        rndx2 = sin(t);
        rndcalcflag = 0;
        return(rndx1 * cos(t));
    }
    else
    {
        rndcalcflag = 1;
        return(rndx2 * rndx1);
    }
}
int rnd(int low, int high)
/* Pick a random integer between low and high */
//int low,high;
{
    int i;
    //float randomperc();

    if(low >= high)
        i = low;
    else
    {
        i = (int)(randomperc() * (high - low + 1)) + low;
        if(i > high) i = high;
    }
    return(i);
}

void warmup_random(float random_seed)
/* Get random off and running */
//float random_seed;
{
    
	int j1, ii;
    double new_random, prev_random;

    oldrand[54] = random_seed;
    new_random = 0.000000001;
    prev_random = random_seed;
    for(j1 = 1 ; j1 <= 54; j1++)
    {
        ii = (21*j1)%54;
        oldrand[ii] = new_random;
        new_random = prev_random-new_random;
        if(new_random<0.0) new_random = new_random + 1.0;
        prev_random = oldrand[ii];
    }

    advance_random();
    advance_random();
    advance_random();

    // v3 (2026): reset the sequence position. randomperc() uses jrand as a
    // rolling index; without the reset, the stream of a fresh warmup_random
    // depends on how many random numbers the PREVIOUS molecule consumed.
    // This was commented out in the original code -- harmless single-threaded
    // (fixed order) but breaks reproducibility when molecules are processed
    // in independent batch chunks (each chunk's first molecule would inherit
    // a different jrand offset).
    jrand = 0;

}

double noise(double mu ,double sigma)
/* normal noise with specified mean & std dev: mu & sigma */
{
    //double randomnormaldeviate();

    return((randomnormaldeviate()*sigma) + mu);
}

float TimeRandomSeed()
{
	srand((unsigned)clock());
	float r = rand()%1001;
	return r/1000.;
}
/*-------------------------------------------------------------*/
