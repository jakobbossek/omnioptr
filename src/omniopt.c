#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "rand.h"

SEXP prepare_objective_vectors_for_r(population *pop, int popsize, int nobj) {
  SEXP ret = PROTECT(allocVector(REALSXP, popsize * nobj));
  double *pret = REAL(ret);
  int k = 0;
  for (int i = 0; i < popsize; ++i) {
    for (int j = 0; j < nobj; ++j) {
      pret[k] = pop->ind[i].obj[j];
      // Rprintf("Obj %i: %.3f ", j, pop->ind[i].obj[j]);
      ++k;
    }
  }
  UNPROTECT(1);
  return ret;
}

SEXP prepare_real_decision_vectors_for_r(population *pop, int popsize, int nreal) {
  SEXP ret = PROTECT(allocVector(REALSXP, popsize * nreal));
  double *pret = REAL(ret);
  int k = 0;
  for (int i = 0; i < popsize; ++i) {
    for (int j = 0; j < nreal; ++j) {
      pret[k] = pop->ind[i].xreal[j];
      // Rprintf("Gene %i: %.3f ", j, pop->ind[i].xreal[j]);
      ++k;
    }
  }
  UNPROTECT(1);
  return ret;
}

SEXP prepare_binary_decision_vectors_for_r(population *pop, int popsize, int nbin) {
  SEXP ret = PROTECT(allocVector(INTSXP, popsize * nbin));
  int *pret = INTEGER(ret);
  int k = 0;
  for (int i = 0; i < popsize; ++i) {
    for (int j = 0; j < nbin; ++j) {
      pret[k] = pop->ind[i].gene[j][0];
      // Rprintf("Gene %i: %i ", j, pop->ind[i].gene[j][0]);
      ++k;
    }
  }
  UNPROTECT(1);
  return ret;
}


SEXP prepare_population_for_r(population *pop, int popsize, int nobj, int nreal, int nbin) {
  SEXP rpop = PROTECT(allocVector(VECSXP, 2));
  if (nreal != 0) {
    SET_VECTOR_ELT(rpop, 0, prepare_real_decision_vectors_for_r(pop, popsize, nreal));
  } else {
    SET_VECTOR_ELT(rpop, 0, prepare_binary_decision_vectors_for_r(pop, popsize, nbin));
  }
  SET_VECTOR_ELT(rpop, 1, prepare_objective_vectors_for_r(pop, popsize, nobj));
  UNPROTECT(1);
  return(rpop);
}

SEXP omnioptC(
  SEXP fun, // smoof function
  SEXP nobjSEXP, // (integer) number of objectives
  SEXP nrealSEXP, // (integer) number of real-valued input dimensions (decision space)
  SEXP min_realvarSEXP, // (double) Vector of lower box constraints
  SEXP max_realvarSEXP, // (double) Vector of lower box constraints
  SEXP nbinSEXP, // (integer) number of binary valued input dimensions
  SEXP nbitsSEXP, // (integervector) number of bits per binary variable
  SEXP pcross_binSEXP, // (double) probability of crossover for binary variables
  SEXP pmut_binSEXP, // (double) probability of mutation for binary variables
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
  SEXP verboseSEXP, // (integer) print messages? (0 or 1)
  SEXP rho) // (environment) used for calling R smoof function from C
{

    popsize = asInteger(popsizeSEXP); // population size (must be multiple of 4)
    ngen = asInteger(ngenSEXP); // number of generations (termination condition)
    nobj = asInteger(nobjSEXP); // number of objectives
    nreal = asInteger(nrealSEXP); // number of real variables
    nbin = asInteger(nbinSEXP); // number of binary variables
    nbits = INTEGER(nbitsSEXP); // number of bits per binary variable
    // for (int i = 0; i < nbin; ++i) {
    //   Rprintf("nBits for %i: %i", i, nbits[i]);
    // }

    //FIXME: Jakob: currently we support only unconstrained optimisation
    ncon = 0; // number of constraints

    // Jakob:  there is no analysis mode from R
    run_mode = 1; // simulation mode, 0 for Analysis mode, 1 for Turbo mode

    // RNG seed
    seed = asReal(seedSEXP); // random number generator seed (real number in [0,1])

    // be verbose?
    int verbose = asInteger(verboseSEXP);

    min_obj = (double *)malloc(nobj*sizeof(double));
    max_obj = (double *)malloc(nobj*sizeof(double));
    epsilon = (double *)malloc(nobj*sizeof(double));

    // box constraints for real variables
    if (nreal > 0) {
      min_realvar = (double *)malloc(nreal*sizeof(double));
      max_realvar = (double *)malloc(nreal*sizeof(double));
      for (int i=0; i<nreal; i++) {
        min_realvar[i] = REAL(min_realvarSEXP)[i];
        max_realvar[i] = REAL(max_realvarSEXP)[i];
      }
    }

    if (nbin > 0) {
      min_binvar = (double *)malloc(nbin*sizeof(double));
      max_binvar = (double *)malloc(nbin*sizeof(double));
      for (int i=0; i<nbin; i++) {
        min_binvar[i] = INF;
        max_binvar[i] = -INF;
      }
    }

    pcross_real = asReal(pcross_realSEXP); // probability of crossover of real variable (0.6-1.0)
    pmut_real = asReal(pmut_realSEXP); // probablity of mutation of real variables (1/nreal)  FIXME: 1/nreal should be default?
    eta_c = asReal(eta_cSEXP); // value of distribution index for crossover (5-20)
    eta_m = asReal(eta_mSEXP); // value of distribution index for mutation (5-50)

    pcross_bin = asReal(pcross_binSEXP); // probability of crossover of binary variable (0.6-1.0)
    pmut_bin = asReal(pmut_binSEXP); // probablity of mutation of binary variables (1/nbin)  FIXME: 1/nreal should be default?

    mate = asInteger(mateSEXP); // choice for selection restriction, 0 for normal selection, 1 for restricted selection
    delta = asReal(deltaSEXP); // delta (0.0 - 1.0) for loose domination

    var_option = asInteger(var_optionSEXP); // variable space niching, 0 for NO, 1 for YES
    obj_option = asInteger(obj_optionSEXP); // objective space niching, 0 for NO, 1 for YES

    input_type = asInteger(input_typeSEXP); // choice for population initialization, 0 for random, 1 for latin-hypercube based sampling, 2 for reading initial population from a file
    frequency = asInteger(frequencySEXP); // frequency with which the population information (only the variables) is to be stored

    char *s;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;


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

    if (verbose == 1) {
      Rprintf("Omni-Optimizer (nreal: %i, nbin: %i, nobj: %i)\n", nreal, nbin, nobj);
      Rprintf("Population size: %i\n", popsize);
      Rprintf("Nr. of generations: %i\n", ngen);
    }

    // <Jakob-start>
    // Many thanks to: https://ro-che.info/articles/2017-08-18-call-r-function-from-c
    // evaluate_pop (parent_pop);

    for (int i = 0; i < popsize; i++)
    {
      if (nreal > 0) {
        double *xreal = parent_pop->ind[i].xreal;

        // convert to SEXP
        SEXP xrealr = PROTECT(allocVector(REALSXP, nreal));
        for (int k = 0; k < nreal; ++k) {
          REAL(xrealr)[k] = xreal[k];
        }

        SEXP call = PROTECT(LCONS(fun, LCONS(xrealr, R_NilValue)));
        SEXP retu = R_forceAndCall(call, 1, rho);

        // update indiviual
        for (int j = 0; j < nobj; ++j) {
          parent_pop->ind[i].obj[j] = REAL(retu)[j];
        }

        // we do not deal with constraint problems
        parent_pop->ind[i].constr_violation = 0.0;

        // drop garbage collector protection
        UNPROTECT(2);
      } else { // nbin > 0
        Rprintf("Evaluating ...");
        // convert to SEXP
        SEXP xbinr = PROTECT(allocVector(INTSXP, nbin));
        for (int k = 0; k < nbin; ++k) {
          INTEGER(xbinr)[k] = parent_pop->ind[i].gene[k][0];
        }

        SEXP call = PROTECT(LCONS(fun, LCONS(xbinr, R_NilValue)));
        SEXP retu = R_forceAndCall(call, 1, rho);

        // update indiviual
        for (int j = 0; j < nobj; ++j) {
          parent_pop->ind[i].obj[j] = REAL(retu)[j];
        }

        // we do not deal with constraint problems
        parent_pop->ind[i].constr_violation = 0.0;

        // drop garbage collector protection
        UNPROTECT(2);
      }
    }
    // <Jakob-end>

    define_epsilon (parent_pop, popsize, epsilon);
    assign_rank_and_crowding_distance (parent_pop);

    int nsaves = ceil(ngen / frequency);
    SEXP r_populations = PROTECT(allocVector(VECSXP, nsaves + 1));
    int curr = 0;

    // save initial population
    SET_VECTOR_ELT(r_populations, curr, prepare_population_for_r(parent_pop, popsize, nobj, nreal, nbin));
    curr++;

    if (ngen==1)
    {
    }
    else
    {
        for (int cgen = 2; cgen<=ngen; cgen++)
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
              if (nreal > 0) {
                double *xreal = child_pop->ind[i].xreal;

                // convert to SEXP
                SEXP xrealr = PROTECT(allocVector(REALSXP, nreal));
                for (int k = 0; k < nreal; ++k) {
                  REAL(xrealr)[k] = xreal[k];
                }

                SEXP call = PROTECT(LCONS(fun, LCONS(xrealr, R_NilValue)));
                SEXP retu = R_forceAndCall(call, 1, rho);

                // update indiviual
                for (int j = 0; j < nobj; ++j) {
                  child_pop->ind[i].obj[j] = REAL(retu)[j];
                }

                // we do not deal with constraint problems
                child_pop->ind[i].constr_violation = 0.0;

                // drop garbage collector protection
                UNPROTECT(2);
              } else { // nbin > 0

                // convert to SEXP
                SEXP xbinr = PROTECT(allocVector(INTSXP, nbin));
                for (int k = 0; k < nbin; ++k) {
                  INTEGER(xbinr)[k] = parent_pop->ind[i].gene[k][0];
                }

                SEXP call = PROTECT(LCONS(fun, LCONS(xbinr, R_NilValue)));
                SEXP retu = R_forceAndCall(call, 1, rho);

                // update indiviual
                for (int j = 0; j < nobj; ++j) {
                  parent_pop->ind[i].obj[j] = REAL(retu)[j];
                }

                // we do not deal with constraint problems
                parent_pop->ind[i].constr_violation = 0.0;

                // drop garbage collector protection
                UNPROTECT(2);
              }
            }
            // <Jakob-end>

            merge (parent_pop, child_pop, mixed_pop);
            define_epsilon (mixed_pop, 2*popsize, epsilon);
            fill_nondominated_sort (mixed_pop, parent_pop);

            if (verbose && (cgen % 10 == 0)) {
              Rprintf("Generation %i finished.\n", cgen);
            }

            if (cgen % frequency == 0 && cgen >= frequency) {
              SET_VECTOR_ELT(r_populations, curr, prepare_population_for_r(parent_pop, popsize, nobj, nreal, nbin));
              curr++;
            }
      }
    }

    if (verbose == 1) {
      Rprintf("Routine successfully terminated.\n");
    }

    // Generate list to return both vectors simultaneously
    SEXP rout = PROTECT(allocVector(VECSXP, 3));
    if (nreal != 0) {
      SET_VECTOR_ELT(rout, 0, prepare_real_decision_vectors_for_r(parent_pop, popsize, nreal));
    } else {
      SET_VECTOR_ELT(rout, 0, prepare_binary_decision_vectors_for_r(parent_pop, popsize, nbin));
    }
    SET_VECTOR_ELT(rout, 1, prepare_objective_vectors_for_r(parent_pop, popsize, nobj));
    SET_VECTOR_ELT(rout, 2, r_populations);

    Rprintf("Fuck");

    free(s);
    if (nreal!=0) {
      free (min_realvar);
      free (max_realvar);
    }
    if (nbits!=0) {
      free (min_binvar);
      free (max_binvar);
    }

    // clean up
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2 * popsize);
    free(parent_pop);
    free(child_pop);
    free(mixed_pop);

    UNPROTECT(2);
    return rout;
}
