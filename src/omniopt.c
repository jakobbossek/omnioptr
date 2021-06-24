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

SEXP omnioptC(
  SEXP fun, // smoof function
  SEXP nobjSEXP, // (integer) number of objectives
  SEXP nrealSEXP, // (integer) number of input dimensions (decision space)
  SEXP popsizeSEXP, // (integer) population size
  SEXP ngenSEXP, // (integer) number of generations
  SEXP pcross_realSEXP, // (double) probability of crossover
  SEXP pmut_realSEXP, // (double) probablity of mutation of real variables (1/nreal)
  SEXP eta_cSEXP, // (double) value of distribution index for crossover (5-20)
  SEXP eta_mSEXP, // (double) value of distribution index for mutation (5-50)
  SEXP mateSEXP, // (integer) choice for selection restriction, 0 for normal selection, 1 for restricted selection
  SEXP deltaSEXP, // (double) delta (0.0 - 1.0) for loose domination
  SEXP var_optionSEXP, // (integer) variable space niching, 0 for NO, 1 for YES
  SEXP obj_optionSEXP, // (integer) objective space niching, 0 for NO, 1 for YES
  SEXP input_typeSEXP, // (integer) choice for population initialization, 0 for random, 1 for latin-hypercube based sampling, [JAKOB]: unsupoported third option: 2 for reading initial population from a file
  SEXP frequencySEXP, // frequency with which the population information (only the variables) is to be stored
  SEXP seedSEXP, // (double) RNG seed in (0,1)
  SEXP rho) // (environment) used for calling R smoof function from C
{

    popsize = asInteger(popsizeSEXP); // population size (must be multiple of 4)
    ngen = asInteger(ngenSEXP); // number of generations (termination condition)
    nobj = asInteger(nobjSEXP); // number of objectives
    nreal = asInteger(nrealSEXP); // number of real variables

    //FIXME: Jakob: currently we support only continuous functions
    ncon = 0; // number of constraints
    nbin = 0; // number of binary variables

    // Jakob:  there is no analysis mode from R
    run_mode = 1; // simulation mode, 0 for Analysis mode, 1 for Turbo mode

    // RNG seed
    seed = asReal(seedSEXP); // random number generator seed (real number in [0,1])

    min_obj = (double *)malloc(nobj*sizeof(double));
    max_obj = (double *)malloc(nobj*sizeof(double));
    epsilon = (double *)malloc(nobj*sizeof(double));

    // box constraints for real variables
    //FIXME: Jakob: hardcoded at the moment
    min_realvar = (double *)malloc(nreal*sizeof(double));
    max_realvar = (double *)malloc(nreal*sizeof(double));
    int i;
    for (i=0; i<nreal; i++) {
      min_realvar[i] = 0;
      max_realvar[i] = 1;
    }

    pcross_real = 0.5; // probability of crossover of real variable (0.6-1.0)
    pmut_real = 0.5; // probablity of mutation of real variables (1/nreal)  FIXME: 1/nreal should be default?
    eta_c = 20; // value of distribution index for crossover (5-20)
    eta_m = 20; // value of distribution index for mutation (5-50)

    mate = 0; // choice for selection restriction, 0 for normal selection, 1 for restricted selection
    delta = 0.001; // delta (0.0 - 1.0) for loose domination

    var_option = 0; // variable space niching, 0 for NO, 1 for YES
    obj_option = 1; // objective space niching, 0 for NO, 1 for YES

    input_type = 0; // choice for population initialization, 0 for random, 1 for latin-hypercube based sampling, 2 for reading initial population from a file
    frequency = 1; // frequency with which the population information (only the variables) is to be stored

    // // Jakob: do all these sanity checks (essentially taken from src/input.c)
    // if (popsize<4 || (popsize%4)!= 0) {
    //     Rprintf("\n population size read is : %d",popsize);
    //     Rprintf("\n Wrong population size entered, hence exiting \n");
    //     exit (1);
    // }
    // if (ngen<1) {
    //     Rprintf("\n number of generations read is : %d",ngen);
    //     Rprintf("\n Wrong nuber of generations entered, hence exiting \n");
    //     exit (1);
    // }
    // if (nobj<1) {
    //     Rprintf("\n number of objectives entered is : %d",nobj);
    //     Rprintf("\n Wrong number of objectives entered, hence exiting \n");
    //     exit (1);
    // }
    // if (ncon<0) {
    //     Rprintf("\n number of constraints entered is : %d",ncon);
    //     Rprintf("\n Wrong number of constraints enetered, hence exiting \n");
    //     exit (1);
    // }
    // if (nreal<0) {
    //     Rprintf("\n number of real variables entered is : %d",nreal);
    //     Rprintf("\n Wrong number of variables entered, hence exiting \n");
    //     exit (1);
    // }
    // if (nreal != 0) {
    //     if (pcross_real<0.0 || pcross_real>1.0)
    //     {
    //         Rprintf("\n Probability of crossover entered is : %e",pcross_real);
    //         Rprintf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
    //         exit (1);
    //     }
    //     if (pmut_real<0.0 || pmut_real>1.0)
    //     {
    //         Rprintf("\n Probability of mutation entered is : %e",pmut_real);
    //         Rprintf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
    //         exit (1);
    //     }
    //     if (eta_c<=0)
    //     {
    //         Rprintf("\n The value entered is : %e",eta_c);
    //         Rprintf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
    //         exit (1);
    //     }
    //     // Rprintf ("\n Enter the value of distribution index for mutation (5-50): ");
    //     // scanf ("%lf",&eta_m);
    //     if (eta_m<=0)
    //     {
    //         Rprintf("\n The value entered is : %e",eta_m);
    //         Rprintf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
    //         exit (1);
    //     }
    // }
    // if (nreal==0 && nbin==0) {
    //     Rprintf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
    //     exit(1);
    // }
    // if (mate!=0 && mate!=1) {
    //     Rprintf("\n Selection restriction option read as : %d",mate);
    //     Rprintf("\n Wrong choice entered for selection restriction, hence exiting \n");
    //     exit(1);
    // }

    // if (delta<0.0 || delta>1.0) {
    //     Rprintf("\n Wrong value entered for delta for loose domination, the value read was : %e",delta);
    //     Rprintf("\n Exiting \n");
    //     exit(1);
    // }
    // if (var_option!=0 && var_option!=1) {
    //     Rprintf("\n The choice read for variable space niching is %d",var_option);
    //     Rprintf("\n Wrong choice entered, hence exiting \n");
    // }
    // if (obj_option!=0 && obj_option!=1) {
    //     Rprintf("\n The choice read for objective space niching is %d",obj_option);
    //     Rprintf("\n Wrong choice entered, hence exiting \n");
    //     exit(1);
    // }
    // if (var_option==0 && obj_option==0) {
    //     Rprintf("\n Both variable and objective space niching cannot be zero, exiting\n");
    //     exit(1);
    // }
    // if (input_type!=0 && input_type!=1 && input_type!=2) {
    //     Rprintf("\n Wrong choice entered for population initialization, choice entered was %d",input_type);
    //     Rprintf("\n Exiting \n");
    //     exit(1);
    // }
    // if (frequency<1 || frequency>ngen) {
    //     Rprintf("\n Wrong value of frequency entered, the value read is : %d",frequency);
    //     Rprintf("\n It should be in the range (1 - pop_size), exiting \n");
    //     exit(1);
    // }
    // if (run_mode!=0 && run_mode!=1) {
    //     Rprintf("\n Value read for simulation mode is : %d",run_mode);
    //     Rprintf("\n Wrong value entered, hence exiting \n");
    //     exit(1);
    // }

    char *s;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

    if (seed<=0.0 || seed>=1.0)
    {
        Rprintf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    pdefinit();
    s = (char *)malloc(100*sizeof(double));
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    randomize();
    initialize_pop (parent_pop);
    decode_pop(parent_pop);
    Rprintf("Initialization done, now performing first generation.\n");

    // <Jakob-start>
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
        Rprintf("O%i: %.3f\n", k, REAL(retu)[k]);
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
    // <Jakob-end>

    define_epsilon (parent_pop, popsize, epsilon);
    assign_rank_and_crowding_distance (parent_pop);
    // report_pop (parent_pop, fpt1);

    if (ngen==1)
    {
    // Jakob: we do not need this
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
            // evaluate_pop(child_pop);

            // <Jakob-start>
            // FIXME: copy and paste from above
            // Many thanks to: https://ro-che.info/articles/2017-08-18-call-r-function-from-c

            for (int i = 0; i < popsize; i++)
            {
              double *xreal = child_pop->ind[i].xreal;

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
                Rprintf("O%i: %.3f\n", k, REAL(retu)[k]);
              }

              // update indiviual
              for (int j = 0; j < nobj; ++j) {
                child_pop->ind[i].obj[j] = REAL(retu)[j];
              }

              // we do not deal with constraint problems
              child_pop->ind[i].constr_violation = 0.0;

              // drop garbage collector protection
              UNPROTECT(2);
            }
            // <Jakob-end>
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
    }
    if (nreal!=0) {
        free (min_realvar);
        free (max_realvar);
    }
    if (nbin!=0) {
        free (nbits);
        free (min_binvar);
        free (max_binvar);
    }
    free(s);

    // Jakob: we need to construct a vector and convert to matrix in R
    //FIXME: I guess parent_pop is what we need here
    SEXP r_pareto_set = PROTECT(allocVector(REALSXP, nreal * popsize));
    SEXP r_pareto_front = PROTECT(allocVector(REALSXP, nobj * popsize));

    // Prepare Pareto-set approximation
    int k = 0;
    for (i = 0; i < popsize; ++i) {
      for (int j = 0; j < nreal; ++j) {
        REAL(r_pareto_set)[k] = parent_pop->ind[i].xreal[j];
        ++k;
      }
    }

    // Prepare Pareto-front approximation
    k = 0;
    for (i = 0; i < popsize; ++i) {
      for (int j = 0; j < nobj; ++j) {
        REAL(r_pareto_front)[k] = parent_pop->ind[i].obj[j];
        ++k;
      }
    }

    // Generate list to return both vectors simultaneously
    SEXP rout = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rout, 0, r_pareto_set);
    SET_VECTOR_ELT(rout, 1, r_pareto_front);

    // clean up
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2 * popsize);
    free(parent_pop);
    free(child_pop);
    free(mixed_pop);

    Rprintf("\n Routine successfully exited \n");
    UNPROTECT(3);
    return rout;
}
