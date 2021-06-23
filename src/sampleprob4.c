/* Test problem definitions */
/* Compilation of various test problems taken from literature */
/* Self generated and from non-standard sources (no references available) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
# include "pdef.h"

/* # define prob1 */
/* # define prob2 */
/* # define prob3 */
/* # define prob4 */
/* # define atrr_13_1 */
/* # define atrr_13_2 */
/* # define atrr_14 */
/* # define atrr_15 */
/* # define atrr_16 */

/*  Test problem PROB1
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef prob1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = pow (sin(xreal[0]*PI), 2.0);
    return;
}
#endif

/*  Test problem PROB2
    # of real variables = 5
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef prob2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
	obj[0] = sin(xreal[0]*PI)+sin(xreal[1]*PI)+sin(xreal[2]*PI)+sin(xreal[3]*PI)+sin(xreal[4]*PI);
	obj[1] = cos(xreal[0]*PI)+cos(xreal[1]*PI)+cos(xreal[2]*PI)+cos(xreal[3]*PI)+cos(xreal[4]*PI);
	return;
}
#endif

/*  Test problem PROB3
    # of real variables = 5
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef prob3
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
	double temp;
	temp = xreal[0]+xreal[1]+xreal[2]+xreal[3]+xreal[4];
	temp = temp*PI;
	obj[0] = sin(temp);
	obj[1] = cos(temp);
	return;
}
#endif

/*  Test problem PROB4
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 0
    */

#ifdef prob4
void test_problem (double *xreal, double *bin, int **gene, double *obj, double *constr)
{
    int i, j;
    double temp;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        temp = 0.0;
        for (j=0; j<i; j++)
        {
            temp += xreal[j];
        }
        obj[0] += square(temp);
    }
    return;
}
#endif

#ifdef atrr_13_1
void test_problem (double *x, double *bb, int **gene, double *obj, double *gg)
{
    double a, b, g;
    a = (1.0-exp(-4.0*x[0]))/(1.0-exp(-4.0));
    b = 2.0 - exp(-2.0*x[1])*cos(8.0*PI*x[1]);
    g = 2.0 - pow((a/b),0.6);
    obj[0] = a;
    obj[1] = g*b;
    return;
}
#endif

#ifdef atrr_13_2
void test_problem (double *x, double *bb, int **gene, double *obj, double *gg)
{
    double a, b, g;
    double x2;
    a = (1.0-exp(-4.0*x[0]))/(1.0-exp(-4.0));
    x2 = 1.0 - 0.1*x[1] - 0.2*square(x[1]) - 0.3*x[2] - 0.1*square(x[2]) - 0.3*x[1]*x[2];
    b = 2.0 - exp(-2.0*x2)*cos(8.0*PI*x2);
    g = 2.0 - pow((a/b),0.6);
    obj[0] = a;
    obj[1] = g*b;
    return;
}
#endif

#ifdef atrr_14
void test_problem (double *x, double *bb, int **gene, double *obj, double *gg)
{
    double a, b, g;
    double x2;
    a = (1.0-exp(-4.0*x[0]))/(1.0-exp(-4.0));
    x2 = 1.0 - 0.1*x[1] - 0.3*x[1]*square(cos(2.0*PI*x[1])) - 0.2*x[2] - 0.2*square(x[2]);
    b = 2.0 - exp(-2.0*x2)*cos(8.0*PI*x2);
    g = 2.0 - pow((a/b),0.6);
    obj[0] = a;
    obj[1] = g*b;
    return;
}
#endif

#ifdef atrr_15
void test_problem (double *x, double *bb, int **gene, double *obj, double *gg)
{
    double a, b, g;
    double x2;
    a = (1.0-exp(-3.0*x[0]))/(1.0-exp(-3.0));
    x2 = 0.2 + 0.2*x[1] + 0.6*square(x[1]);
    b = 3.0 - exp(-x2)*cos(2.0*PI*x2) -exp(-x[1])*cos(4.0*PI*x[1]) ;
    g = 2.0 - pow((a/b),0.4) - (a/b)*cos(8.0*PI*a*a);
    obj[0] = a;
    obj[1] = g*b;
    return;
}
#endif

#ifdef atrr_16
void test_problem (double *x, double *bb, int **gene, double *obj, double *gg)
{
    double a, b, g;
    double x2;
    a = (1.0-exp(-3.0*x[0]))/(1.0-exp(-3.0));
    x2 = 0.2 + 3.396*x[1]-3.396*square(x[1]) - exp(0.0957*x[1])*exp(-square((x[1]-0.5)/0.004));
    b = 3.0 - exp(-x2)*cos(2.0*PI*x2) -exp(-x[1])*cos(4.0*PI*x[1]) ;
    g = 2.0 - pow((a/b),0.4) - (a/b)*cos(8.0*PI*a*a);
    obj[0] = a;
    obj[1] = g*b;
    return;
}
#endif
