/* Test problem definitions */
/* Compilation of various test problems taken from literature */
/* Standard single-objective test problems */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"
#include "pdef.h"

/* # define ackley */
/* # define griewank */
/* # define himmel */
/* # define rana */
/* # define rastrigin */
/* # define rosenbrock */
/* # define schwefel */
/* # define sphere */
/* # define weierstrass */

/*  Test problem ackley
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

# ifdef ackley
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0; i<nreal; i++)
    {
        sum1 += square(x[i]);
        sum2 += cos(2.0*PI*x[i]);
    }
    sum1 = -0.2*sqrt(sum1/nreal);
    sum2 = sum2/nreal;
    obj[0] = 20.0 + exp(1.0) - 20.0*exp(sum1) - exp(sum2);
    return;
}
# endif

/*  Test problem griewank
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

# ifdef griewank
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    double s;
    double p;
    s = 0.0;
    p = 1.0;
    for (i=0; i<nreal; i++)
    {
        s += square(x[i]);
        p = p*cos(x[i]/sqrt((double)(i+1)));
    }
    obj[0] = 1.0 + s/4000.0 - p;
    return;
}
# endif

/*  Test problem himmel
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef himmel
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = pow((xreal[0]*xreal[0] + xreal[1]-11.0),2.0) + pow((xreal[0] + xreal[1]*xreal[1]-7.0),2.0);
    return;
}
#endif

/*  Test problem rana
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

# ifdef rana
void test_problem (double *xreal, double *b, int **gene, double *obj, double *constr)
{
    int i;
    double x, y;
    obj[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        x = xreal[i];
        y = xreal[i+1];
        obj[0] += x*sin(sqrt(fabs(y+1.0-x)))*cos(sqrt(fabs(y+1.0+x))) + (y+1)*sin(sqrt(fabs(x+y+1.0)))*cos(sqrt(fabs(y+1.0-x)));
    }
    return;
}
#endif

/*  Test problem rastrigin
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef rastrigin
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 10.0*nreal;
    for (i=0; i<nreal; i++)
    {
        obj[0] += (xreal[i]*xreal[i] - 10.0*cos(2.0*PI*xreal[i]));
    }
    return;
}
#endif

/*  Test problem rosenbrock
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef rosenbrock
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal-1; i++)
    {
        obj[0] += 100.0*pow((xreal[i]*xreal[i]-xreal[i+1]),2.0) + 1.0*pow((xreal[i]-1.0),2.0);
    }
    return;
}
#endif

/*  Test problem schwefel
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef schwefel
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        obj[0] += -xreal[i]*sin(sqrt(fabs(xreal[i])));
    }
    obj[0] += 418.9829*nreal;
    return;
}
#endif

/*  Test problem sphere
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef sphere
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        obj[0] += xreal[i]*xreal[i];
    }
    return;
}
#endif

/*  Test problem weierstrass
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef weierstrass
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i, k, k_max;
    double a, b;
    a = 0.5;
    b = 3.0;
    k_max = 20;
    obj[0] = 2.0*nreal;
    for (i=0; i<nreal; i++)
    {
        for (k=0; k<=k_max; k++)
        {
            obj[0] += pow(a,k)*cos(2.0*PI*pow(b,k)*xreal[i]);
        }
    }
    return;
}
#endif
