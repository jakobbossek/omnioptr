/* Test problem definitions */
/* Compilation of various test problems taken from literature */
/* Problems taken from Shinn et. al. (many of them are same as given in problemdef2.c) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
# include "pdef.h"

/* # define f1 */
/* # define f2 */
/* # define f3 */
/* # define f4 */
/* # define f5 */
/* # define f6 */
/* # define f7 */
/* # define f8 */
/* # define f9 */
/* # define f10 */
/* # define f11 */
/* # define f12 */

# ifdef f1
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        obj[0] += (sin(x[i]) + sin(2.0*x[i]/3.0));
    }
    return;
}
# endif

# ifdef f2
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<(nreal-1); i++)
    {
        obj[0] += ( sin(x[i]+x[i+1]) + sin(2.0*(x[i]*x[i+1])/3.0));
    }
    return;
}
# endif

# ifdef f3
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    int temp;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        if ((x[i]+0.5)<0.0)
        {
            temp = (x[i]+0.5);
            temp = temp-1;
        }
        else
        {
            temp = (x[i]+0.5);
        }
        obj[0] += (double)square(temp);
    }
    return;
}
# endif

#ifdef f4
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

# ifdef f5
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        obj[0] += square(x[i]);
    }
    return;
}
# endif

# ifdef f6
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        obj[0] = obj[0] - (x[i]*sin(10.0*PI*x[i]));
    }
    return;
}
# endif

# ifdef f7
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    obj[0] = 0.0;
    for (i=0; i<nreal; i++)
    {
        if (fabs(x[i]) < 0.000001)
        {
            obj[0] += 1.0;
        }
        else
        {
            obj[0] += fabs( (sin(10.0*PI*x[i])) / (10.0*PI*x[i]) );
        }
    }
    return;
}
# endif

# ifdef f8
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

#ifdef f9
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

#ifdef f10
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

# ifdef f11
void test_problem (double *x, double *b, int **gene, double *obj, double *constr)
{
    int i;
    int temp;
    obj[0] = 6.0*nreal;
    for (i=0; i<nreal; i++)
    {
        if (x[i] < 0.0)
        {
            temp = x[i];
            temp = temp-1;
        }
        else
        {
            temp = x[i];
        }
        obj[0] += (double)temp;
    }
    return;
}
# endif

# ifdef f12
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
