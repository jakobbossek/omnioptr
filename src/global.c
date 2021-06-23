#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"

int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
double *min_obj;
double *max_obj;
double *epsilon;
int bitlength;
double delta;
int mate;
int input_type;
int run_mode;
int var_option;
int obj_option;
int frequency;
int var1;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;
int choice;