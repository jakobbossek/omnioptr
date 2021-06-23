#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "rand.h"

#include "macros.h"

SEXP R_apply_fun(SEXP f, SEXP x, SEXP rho) {
  SEXP call = PROTECT(LCONS(f, LCONS(x, R_NilValue)));
  SEXP val = R_forceAndCall(call, 1, rho);
  UNPROTECT(1);
  return val;
}

SEXP omnioptC(SEXP fun, SEXP popsizeSEXP, SEXP ngenSEXP, SEXP rho) {

    // double tt[2] = {2, 3.1};
    // int size = 2;
    // SEXP ttr = PROTECT(allocVector(REALSXP, size));
    // for (int i = 0; i < 2; ++i) {
    //   REAL(ttr)[i] = tt[i];
    // }
    // Rprintf( "Program starting...\n" );
    // SEXP call = PROTECT(LCONS(fun, LCONS(ttr, R_NilValue)));
    // Rprintf( "Created call...\n" );
    // SEXP val = R_forceAndCall(call, 1, rho);
    // Rprintf( "Program running...\n" );
    // UNPROTECT(2);
    // return val;
    // return ScalarReal(10);


    int i;


    popsize = asInteger(popsizeSEXP); // population size (must be multiple of 4)
    ngen = asInteger(ngenSEXP); // number of generations (termination condition)
    nobj = 2; // number of objectives
    ncon = 0; // number of constraints
    nreal = 2; // number of real variables
    nbin = 0; // number of binary variables

    min_obj = (double *)malloc(nobj*sizeof(double));
    max_obj = (double *)malloc(nobj*sizeof(double));
    epsilon = (double *)malloc(nobj*sizeof(double));

    // box constraints for real variables
    min_realvar = (double *)malloc(nreal*sizeof(double));
    max_realvar = (double *)malloc(nreal*sizeof(double));
    for (i=0; i<nreal; i++) {
      min_realvar[i] = -5;
      max_realvar[i] = 5;
    }


    pcross_real = 0.5; // probability of crossover of real variable (0.6-1.0)
    pmut_real = 0.5; // probablity of mutation of real variables (1/nreal)  FIXME: 1/nreal should be default?
    eta_c = 10; // value of distribution index for crossover (5-20)
    eta_m = 10; // value of distribution index for mutation (5-50)

    mate = 1; // choice for selection restriction, 0 for normal selection, 1 for restricted selection
    delta = 0.5; // delta (0.0 - 1.0) for loose domination

    var_option = 0; // variable space niching, 0 for NO, 1 for YES
    obj_option = 1; // objective space niching, 0 for NO, 1 for YES

    input_type = 0; // choice for population initialization, 0 for random, 1 for latin-hypercube based sampling, 2 for reading initial population from a file
    frequency = 1; // frequency with which the population information (only the variables) is to be stored
    run_mode = 1; // simulation mode, 0 for Analysis mode, 1 for Turbo mode

    if (popsize<4 || (popsize%4)!= 0)
    {
        printf("\n population size read is : %d",popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    // printf("\n Enter the number of generations : ");
    // scanf("%d",&ngen);
    if (ngen<1)
    {
        printf("\n number of generations read is : %d",ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    // printf("\n Enter the number of objectives : ");
    // scanf("%d",&nobj);
    if (nobj<1)
    {
        printf("\n number of objectives entered is : %d",nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }
    // printf("\n Enter the number of constraints : ");
    // scanf("%d",&ncon);
    if (ncon<0)
    {
        printf("\n number of constraints entered is : %d",ncon);
        printf("\n Wrong number of constraints enetered, hence exiting \n");
        exit (1);
    }
    // printf("\n Enter the number of real variables : ");
    // scanf("%d",&nreal);
    if (nreal<0)
    {
        printf("\n number of real variables entered is : %d",nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }
    if (nreal != 0)
    {
        // min_realvar = (double *)malloc(nreal*sizeof(double));
        // max_realvar = (double *)malloc(nreal*sizeof(double));
        // for (i=0; i<nreal; i++)
        // {
        //     printf ("\n Enter the lower limit of real variable %d : ",i+1);
        //     scanf ("%lf",&min_realvar[i]);
        //     printf ("\n Enter the upper limit of real variable %d : ",i+1);
        //     scanf ("%lf",&max_realvar[i]);
        //     if (max_realvar[i] <= min_realvar[i])
        //     {
        //         printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
        //         exit(1);
        //     }
        // }
        // printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        // scanf ("%lf",&pcross_real);
        if (pcross_real<0.0 || pcross_real>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_real);
            printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        // printf ("\n Enter the probablity of mutation of real variables (1/nreal) : ");
        // scanf ("%lf",&pmut_real);
        if (pmut_real<0.0 || pmut_real>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        // printf ("\n Enter the value of distribution index for crossover (5-20): ");
        // scanf ("%lf",&eta_c);
        if (eta_c<=0)
        {
            printf("\n The value entered is : %e",eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }
        // printf ("\n Enter the value of distribution index for mutation (5-50): ");
        // scanf ("%lf",&eta_m);
        if (eta_m<=0)
        {
            printf("\n The value entered is : %e",eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }
    }
    // printf("\n Enter the number of binary variables : ");
    // scanf("%d",&nbin);
    // if (nbin<0)
    // {
    //     printf ("\n number of binary variables entered is : %d",nbin);
    //     printf ("\n Wrong number of binary variables entered, hence exiting \n");
    //     exit(1);
    // }
    // if (nbin != 0)
    // {
    //     nbits = (int *)malloc(nbin*sizeof(int));
    //     min_binvar = (double *)malloc(nbin*sizeof(double));
    //     max_binvar = (double *)malloc(nbin*sizeof(double));
    //     for (i=0; i<nbin; i++)
    //     {
    //         printf ("\n Enter the number of bits for binary variable %d : ",i+1);
    //         scanf ("%d",&nbits[i]);
    //         if (nbits[i] < 1)
    //         {
    //             printf("\n Wrong number of bits for binary variable entered, hence exiting");
    //             exit(1);
    //         }
    //         printf ("\n Enter the lower limit of binary variable %d : ",i+1);
    //         scanf ("%lf",&min_binvar[i]);
    //         printf ("\n Enter the upper limit of binary variable %d : ",i+1);
    //         scanf ("%lf",&max_binvar[i]);
    //         if (max_binvar[i] <= min_binvar[i])
    //         {
    //             printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
    //             exit(1);
    //         }
    //     }
    //     printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
    //     scanf ("%lf",&pcross_bin);
    //     if (pcross_bin<0.0 || pcross_bin>1.0)
    //     {
    //         printf("\n Probability of crossover entered is : %e",pcross_bin);
    //         printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
    //         exit (1);
    //     }
    //     printf ("\n Enter the probability of mutation of binary variables (1/nbits): ");
    //     scanf ("%lf",&pmut_bin);
    //     if (pmut_bin<0.0 || pmut_bin>1.0)
    //     {
    //         printf("\n Probability of mutation entered is : %e",pmut_bin);
    //         printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
    //         exit (1);
    //     }
    // }
    if (nreal==0 && nbin==0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }
    // printf("\n Enter the choice for selection restriction, 0 for normal selection, 1 for restricted selection : ");
    // scanf("%d",&mate);
    if (mate!=0 && mate!=1)
    {
        printf("\n Selection restriction option read as : %d",mate);
        printf("\n Wrong choice entered for selection restriction, hence exiting \n");
        exit(1);
    }
    // printf("\n Enter the value of delta (0.0 - 1.0) for loose domination : ");
    // scanf("%lf",&delta);
    if (delta<0.0 || delta>1.0)
    {
        printf("\n Wrong value entered for delta for loose domination, the value read was : %e",delta);
        printf("\n Exiting \n");
        exit(1);
    }
    // printf("\n Do you want to enable variable space niching, 0 for NO, 1 for YES : ");
    // scanf("%d",&var_option);
    if (var_option!=0 && var_option!=1)
    {
        printf("\n The choice read for variable space niching is %d",var_option);
        printf("\n Wrong choice entered, hence exiting \n");
        exit(1);
    }
    // printf("\n Do you want to enable objective space niching, 0 for NO, 1 for YES : ");
    // scanf("%d",&obj_option);
    if (obj_option!=0 && obj_option!=1)
    {
        printf("\n The choice read for objective space niching is %d",obj_option);
        printf("\n Wrong choice entered, hence exiting \n");
        exit(1);
    }
    if (var_option==0 && obj_option==0)
    {
        printf("\n Both variable and objective space niching cannot be zero, exiting\n");
        exit(1);
    }
    // printf("\n Enter the choice for population initialization, 0 for random, 1 for latin-hypercube based sampling, 2 for reading initial population from a file : ");
    // scanf("%d",&input_type);
    if (input_type!=0 && input_type!=1 && input_type!=2)
    {
        printf("\n Wrong choice entered for population initialization, choice entered was %d",input_type);
        printf("\n Exiting \n");
        exit(1);
    }
    // printf("\n Enter the frequency with which the population information (only the variables) is to be stored : ");
    // scanf("%d",&frequency);
    if (frequency<1 || frequency>ngen)
    {
        printf("\n Wrong value of frequency entered, the value read is : %d",frequency);
        printf("\n It should be in the range (1 - pop_size), exiting \n");
        exit(1);
    }
    // printf("\n Enter the simulation mode, 0 for Analysis mode, 1 for Turbo mode : ");
    // scanf("%d",&run_mode);
    if (run_mode!=0 && run_mode!=1)
    {
        printf("\n Value read for simulation mode is : %d",run_mode);
        printf("\n Wrong value entered, hence exiting \n");
        exit(1);
    }


    // double te = 1.0;
    // return ScalarReal(te);
    // FILE *fpt1;
    // FILE *fpt2;
    // FILE *fpt3;
    // FILE *fpt4;
    // FILE *fpt5;
    // FILE *fpt6;
    // FILE *fpt;
    // FILE *fpt7;
    // FILE *gp;
    char *s;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;
    printf("\n Omni-optimizer code developed at KanGAL, IIT Kanpur.");
    printf("\n Copyright (C) 2006 by Dr. Kalyanmoy Deb and Santosh Tiwari");
    printf("\n");
    fflush(stdout);

    seed = 0.5;
    if (seed<=0.0 || seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    // input();
    pdefinit();
    // fpt1 = fopen("initial_pop.out","w");
    // fpt2 = fopen("final_pop.out","w");
    // fpt3 = fopen("nondominated_pop.out","w");
    // fpt4 = fopen("output.out","w");
    // fpt5 = fopen("problem_data.out","w");
    // fpt6 = fopen("plot.out","w");
    // fpt7 = fopen("run.out","a");
    s = (char *)malloc(100*sizeof(double));
    // output (fpt1, fpt2, fpt3, fpt4, fpt5);
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    randomize();
    initialize_pop (parent_pop);
    decode_pop(parent_pop);
    printf("\n Initialization done, now performing first generation");
    fflush(stdout);

    // <JB-start>
    // Many thanks to: https://ro-che.info/articles/2017-08-18-call-r-function-from-c
    evaluate_pop (parent_pop);

    for (int i = 0; i < popsize; i++)
    {
      double *xreal = parent_pop->ind[i].xreal;

      // convert to SEXP
      SEXP xrealr = PROTECT(allocVector(REALSXP, nreal));
      for (int k = 0; k < nreal; ++k) {
        REAL(xrealr)[k] = xreal[k];
      }

      // debug
      for (int k = 0; k < nreal; ++k) {
        Rprintf("%.2f\n", REAL(xrealr)[k]);
      }

      SEXP call = PROTECT(LCONS(fun, LCONS(xrealr, R_NilValue)));
      SEXP retu = R_forceAndCall(call, 1, rho);

      // debug
      for (int k = 0;  k < nobj; ++k) {
        Rprintf("o%i: %.3f\n", k, REAL(retu)[k]);
      }

      // update indiviual
      for (int j = 0; j < nobj; ++j) {
        parent_pop->ind[i].obj[j] = REAL(retu)[j];
      }

      // we do not deal with constraint problems
      parent_pop->ind[i].constr_violation = 0.0;

      // drop garbage collector protection
      UNPROTECT(2);
    }
    // <JB-end>


    double tt[2] = {2, 3.1};
    int size = 2;
    SEXP ttr = PROTECT(allocVector(REALSXP, size));
    for (int i = 0; i < 2; ++i) {
      REAL(ttr)[i] = tt[i];
    }
    Rprintf( "Program starting...\n" );
    SEXP call = PROTECT(LCONS(fun, LCONS(ttr, R_NilValue)));
    Rprintf( "Created call...\n" );
    SEXP val = R_forceAndCall(call, 1, rho);
    Rprintf( "Program running...\n" );
    UNPROTECT(2);
    return val;


    define_epsilon (parent_pop, popsize, epsilon);
    assign_rank_and_crowding_distance (parent_pop);
    // report_pop (parent_pop, fpt1);

  printf("\n gen = 1");
  fflush(stdout);
    // fflush(fpt1);
    // fflush(fpt2);
    // fflush(fpt3);
    // fflush(fpt4);
    // fflush(fpt5);
    if (ngen==1)
    {
        printf("\n Generations finished, now reporting solutions");
        // report_pop(parent_pop,fpt2);
        // report_feasible(parent_pop,fpt3);
        // report_obj(parent_pop,fpt6);
        // fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
        // fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
        // fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
        // fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
    else
    {
        for (i=2; i<=ngen; i++)
        {
            if (mate==0)
            {
                selection (parent_pop, child_pop);
            }
            else
            {
                restricted_selection(parent_pop, child_pop);
            }
            mutation_pop (child_pop);
            decode_pop(child_pop);
            evaluate_pop(child_pop);
            merge (parent_pop, child_pop, mixed_pop);
            define_epsilon (mixed_pop, 2*popsize, epsilon);
            fill_nondominated_sort (mixed_pop, parent_pop);

            // if (i%frequency==0 && i>=frequency)
            // {
            //     sprintf(s,"pop_var_%d_.out",i);
            //     fpt = fopen(s,"w");
            //     report_var(parent_pop,fpt);
            //     fflush(fpt);
            //     fclose(fpt);
            // }
      printf("\n gen = %d",i);
        }
        printf("\n Generations finished, now reporting solutions");
      //   report_pop(parent_pop,fpt2);
      //   report_feasible(parent_pop,fpt3);
        // report_obj(parent_pop,fpt6);
      //   fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
      //   fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
      //   fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
      //   fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
    // fprintf(fpt7,"%f\t%e\n",seed,parent_pop->ind[0].obj[0]);
    // fflush(fpt7);
    // fclose(fpt7);
    // fflush(stdout);
    // fflush(fpt1);
    // fflush(fpt2);
    // fflush(fpt3);
    // fflush(fpt4);
    // fflush(fpt5);
    // fflush(fpt6);
    // fclose(fpt1);
    // fclose(fpt2);
    // fclose(fpt3);
    // fclose(fpt4);
    // fclose(fpt5);
    // fclose(fpt6);
    // free (epsilon);
    // free (min_obj);
    // free (max_obj);
    if (nreal!=0)
    {
        free (min_realvar);
        free (max_realvar);
    }
    if (nbin!=0)
    {
        free (nbits);
        free (min_binvar);
        free (max_binvar);
    }
    free(s);

    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    free(parent_pop);
    free(child_pop);
    free(mixed_pop);
    printf("\n Routine successfully exited \n");
    return ScalarReal(10.0);
}
