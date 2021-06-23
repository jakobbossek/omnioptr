#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"

double seed;
double oldrand[55];
int jrand;

void randomize()
{
      int j1;
      for(j1=0; j1<=54; j1++)
      {
            oldrand[j1] = 0.0;
      }
      jrand=0;
      warmup_random (seed);
      return;
}

void warmup_random (double seed)
{
      int j1, ii;
      double new_random, prev_random;
      oldrand[54] = seed;
      new_random = 0.000000001;
      prev_random = seed;
      for(j1=1; j1<=54; j1++)
      {
            ii = (21*j1)%54;
            oldrand[ii] = new_random;
            new_random = prev_random-new_random;
            if(new_random<0.0)
            {
                  new_random += 1.0;
            }
            prev_random = oldrand[ii];
      }
      advance_random ();
      advance_random ();
      advance_random ();
      jrand = 0;
      return;
}

void advance_random ()
{
      int j1;
      double new_random;
      for(j1=0; j1<24; j1++)
      {
            new_random = oldrand[j1]-oldrand[j1+31];
            if(new_random<0.0)
            {
                  new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
      }
      for(j1=24; j1<55; j1++)
      {
            new_random = oldrand[j1]-oldrand[j1-24];
            if(new_random<0.0)
            {
                  new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
      }
}

double randomperc()
{
      jrand++;
      if(jrand>=55)
      {
            jrand = 1;
            advance_random();
      }
      return((double)oldrand[jrand]);
}

int rnd (int low, int high)
{
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        res = low + (int)(randomperc()*(high-low+1));
        if (res > high)
        {
            res = high;
        }
    }
    return (res);
}

double rndreal (double low, double high)
{
    return (low + (high-low)*randomperc());
}
