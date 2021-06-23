/* Test problem definitions */
/* Compilation of various test problems taken from literature */
/* Some constrained test problems taken from Dr. Deb's constrainted test problem paper*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
# include "pdef.h"

/* # define test1 */
/* # define test2 */
/* # define test3 */
/* # define test4 */
/* # define test5 */
/* # define test6 */
/* # define test7 */
/* # define test8 */
/* # define test9 */
/* # define test10 */
/* # define test11 */
/* # define test12 */
/* # define beam */

/*  Test problem test1
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 2
    */

#ifdef test1
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = square(x[0]*x[0] + x[1] - 11.0) + square(x[0] + x[1]*x[1] - 7.0);
    g[0] = (square(x[0]) + square(x[1] - 2.5) - 4.84)/4.84;
    g[1] = (4.84 - square(x[0] - 0.05) - square(x[1] - 2.5))/4.84;
    return;
}
#endif

/*  Test problem test2
    # of real variables = 13
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 9
    */

#ifdef test2
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    int i;
    obj[0] = 5.0*(x[0] + x[1] + x[2] + x[3]);
    for (i=0; i<4; i++)
    {
        obj[0] = obj[0] - 5.0*x[i]*x[i];
    }
    for (i=4; i<13; i++)
    {
        obj[0] = obj[0] - x[i];
    }
    g[0] = (10.0 - 2.0*(x[0] + x[1]) - x[9] - x[10] )/10.0;
    g[1] = (10.0 - 2.0*(x[0] + x[2]) - x[9] - x[11] )/10.0;
    g[2] = (10.0 - 2.0*(x[1] + x[2]) - x[10] - x[11] )/10.0;
    g[3] = x[0]  - x[9]/8.0;
    g[4] = x[1]  - x[10]/8.0;
    g[5] = x[2]  - x[11]/8.0;
    g[6] = x[3] + 0.5*x[4] - 0.5*x[9];
    g[7] = x[5] + 0.5*x[6] - 0.5*x[10];
    g[8] = x[7] + 0.5*x[8] - 0.5*x[11];
    return;
}
#endif

/*  Test problem test3
    # of real variables = 8
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 6
    */

#ifdef test3
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = x[0] + x[1] + x[2];
    g[0] = 1.0 - 0.0025*(x[3] + x[5]);
    g[1] = 1.0 - 0.0025*(x[4] + x[6] - x[3]);
    g[2] = 1.0 - 0.01*(x[7] - x[4]);
    g[3] = (x[0]*x[5] - 833.33252*x[3] - 100.0*x[0])/83333.333 + 1.0;
    g[4] = x[1]*x[6]/1250.0 - x[4] - x[1]*x[3]/1250.0 + x[3];
    g[5] = (x[2]*x[7] -x[2]*x[4] + 2500.0*x[4])/1250000.0 - 1.0;
    return;
}
#endif

/*  Test problem test4
    # of real variables = 7
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 4
    */

#ifdef test4
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = square(x[0] - 10.0) + 5.0*square(x[1] -12.0) + square(square(x[2]));
    obj[0] = obj[0] + 3.0*square(x[3] - 11.0) + 10.0*pow(x[4],6.0) + 7.0*square(x[5]);
    obj[0] = obj[0] + square(square(x[6])) -  4.0*x[5]*x[6] - 10.0*x[5] - 8.0*x[6];
    g[0]= (127.0 - 2.0*square(x[0])-3.0*square(square(x[1])) - x[2] - 4.0*square(x[3]) -5.0*x[4])/127.0;
    g[1]= (282.0 - 7.0*x[0] - 3.0*x[1] - 10.0*square(x[2]) -x[3] +x[4])/282.0;
    g[2] = (196.0 - 23.0*x[0] - square(x[1]) - 6.0*square(x[5]) + 8.0*x[6])/196.0;
    g[3] = -4.0*square(x[0]) -square(x[1]) + 3.0*x[0]*x[1] - 2.0*square(x[2]) - 5.0*x[5];
    g[3] = g[3] + 11.0*x[6];
    g[3] = g[3]/11.0;
    return;
}
#endif

/*  Test problem test5
    # of real variables = 5
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 3
    */

#ifdef test5
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = exp(x[0]*x[1]*x[2]*x[3]*x[4]);
    g[0]= EPS - fabs(square(x[0]) + square(x[1]) + square(x[2]) + square(x[3]) + square(x[4]) - 10.0);
    g[1]= EPS - fabs(x[1]*x[2] - 5.0*x[3]*x[4]);
    g[2]= EPS - fabs(square(x[0])*x[0] + square(x[1])*x[1] + 1.0);
    return;
}
#endif

/*  Test problem test6
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 8
    */

#ifdef test6
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = square(x[0]) + square(x[1]) + x[0]*x[1]  - 14.0*x[0] - 16.0*x[1];
    obj[0] = obj[0] + square(x[2] - 10.0) + 4.0*square(x[3] - 5.0) + square(x[4] - 3.0) + 2.0*square(x[5] - 1.0);
    obj[0] = obj[0] + 5.0*square(x[6]) + 7.0*square(x[7] - 11.0) + 2.0*square(x[8] -10.0) + square(x[9] -7.0) + 45.0;
    g[0] =  (105.0 - 4.0*x[0] - 5.0*x[1] + 3.0*x[6]  - 9.0*x[7])/105.0;
    g[1] =  (-10.0*x[0] + 8.0*x[1] + 17.0*x[6] - 2.0*x[7])/10.0;
    g[2] =  (8.0*x[0] - 2.0*x[1] - 5.0*x[8] + 2.0*x[9] + 12.0)/12.0;
    g[3] =  (120.0 - 3.0*square(x[0] -2.0) - 4.0*square(x[1] -3.0) - 2.0*square(x[2]) + 7.0*x[3] )/120.0;
    g[4] =  (40.0 - 5.0*square(x[0]) - 8.0*x[1] - square(x[2] - 6.0) + 2.0*x[3] )/40.0;
    g[5] =  (6.0*x[5] - square(x[0]) - 2.0*square(x[1] - 2.0) + 2.0*x[0]*x[1] - 14.0*x[4] )/6.0;
    g[6] =  (30.0 -0.5*square(x[0] -8.0) - 2.0*square(x[1] - 4.0) - 3.0*square(x[4] )+ x[5] )/30.0;
    g[7] =  (3.0*x[0] - 6.0*x[1] - 12.0*square(x[8] - 8.0) + 7.0*x[9])/12.0;
    return;
}
#endif

/*  Test problem test7
    # of real variables = 5
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 38
    */

#ifdef test7
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    double y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17;
    double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17;
    double ahock[18] = {0.0, 0.0, 17.505, 11.275, 214.228, 7.458, 0.961, 1.612, 0.146, 107.99, 922.693, 926.832, 18.766, 1072.163, 8961.448, 0.063, 71084.33, 2802713.0};
    double bhock[18] = {0.0, 0.0, 1053.6667, 35.03, 665.585, 584.463, 265.916, 7.046, 0.222, 273.366, 1286.105, 1444.046, 537.141, 3247.039, 26844.086, 0.386, 140000.0, 12146108.0};
    y1 = x[1] + x[2] + 41.6;
    c1 = 0.024 * x[3] - 4.62;
    y2 = 12.5/c1 + 12.0;
    c2 = 0.0003535 * x[0]*x[0] + 0.5311*x[0] + 0.08705 * y2 * x[0];
    c3 = 0.052*x[0] + 78.0 + 0.002377 * y2 * x[0];
    y3 = c2/c3;
    y4 = 19.0*y3;
    c4 = 0.04782*(x[0]-y3) + 0.1956*square(x[0]-y3)/x[1] + 0.6376*y4 + 1.594*y3;
    c5 = 100.0*x[1];
    c6 = x[0] - y3 - y4;
    c7 = 0.95 - c4/c5;
    y5 = c6 * c7;
    y6 = x[0] - y5 - y4 - y3;
    c8 = 0.995 * (y4+y5);
    y7 = c8/y1;
    y8 = c8/3798.0;
    c9 = y7 - 0.0663*y7/y8 - 0.3153;
    y9 = 96.82/c9 +  0.321 * y1;
    y10 = 1.29 * y5 + 1.258 * y4 + 2.29 * y3 + 1.71 * y6;
    y11 = 1.71 * x[0] - 0.452 * y4 + 0.58 * y3;
    c10 = 12.3/752.3;
    c11 = 1.75*y2 * 0.995 * x[0];
    c12 = 0.995 * y10 + 1998.0;
    y12 = c10 * x[0] + c11/c12;
    y13 = c12 - 1.75 * y2;
    y14 = 3623.0 + 64.4 * x[1] + 58.4 * x[2] + 146312.0/(y9+x[4]);
    c13 = 0.995 * y10 + 60.8 * x[1] + 48.0*x[3] -0.1121*y14-5095.0;
    y15 = y13/c13;
    y16 = 148000.0 - 331000.0 * y15  + 40*y13 - 61.0*y15*y13;
    c14 = 2324.0 * y10 - 28740000.0 * y2;
    y17 = 14130000.0 - 1328.0 * y10 - 531.0 * y11 + c14/c12;
    c15 = y13/y15 - y13/0.52;
    c16 = 1.104 - 0.72 * y15;
    c17 = y9 + x[4];
    obj[0] = -5.843e-07*y17 + 1.17e-04*y14 + 2.358e-05*y13  + 1.502e-06*y16  + 0.0321*y12 + 0.004324*y5 + 1.0e-04*c15/c16 + 37.48*y2/c12 + 0.1365;
    g[0] = 1.5*x[1] - x[2];
    g[1] = y1/213.1 - 1.0;
    g[2] = 1.0 - y1/405.23;
    g[3] = y2/ahock[2] - 1.0;
    g[4] = y3/ahock[3] - 1.0;
    g[5] = y4/ahock[4] - 1.0;
    g[6] = y5/ahock[5] - 1.0;
    g[7] = y6/ahock[6] - 1.0;
    g[8] = y7/ahock[7] - 1.0;
    g[9] = y8/ahock[8] - 1.0;
    g[10] = y9/ahock[9] - 1.0;
    g[11] = y10/ahock[10] - 1.0;
    g[12] = y11/ahock[11] - 1.0;
    g[13] = y12/ahock[12] - 1.0;
    g[14] = y13/ahock[13] - 1.0;
    g[15] = y14/ahock[14] - 1.0;
    g[16] = y15/ahock[15] - 1.0;
    g[17] = y16/ahock[16] - 1.0;
    g[18] = y17/ahock[17] - 1.0;
    g[19] = 1.0 - y2/bhock[2];
    g[20] = 1.0 - y3/bhock[3];
    g[21] = 1.0 - y4/bhock[4];
    g[22] = 1.0 - y5/bhock[5];
    g[23] = 1.0 - y6/bhock[6];
    g[24] = 1.0 - y7/bhock[7];
    g[25] = 1.0 - y8/bhock[8];
    g[26] = 1.0 - y9/bhock[9];
    g[27] = 1.0 - y10/bhock[10];
    g[28] = 1.0 - y11/bhock[11];
    g[29] = 1.0 - y12/bhock[12];
    g[30] = 1.0 - y13/bhock[13];
    g[31] = 1.0 - y14/bhock[14];
    g[32] = 1.0 - y15/bhock[15];
    g[33] = 1.0 - y16/bhock[16];
    g[34] = 1.0 - y17/bhock[17];
    g[35] = y4 - 0.28/(0.72*y5);
    g[36] = 1.0 - 3496.0 * y2/(21.0 * c12);
    g[37] = (62212.0/c17 - y1)/110.6 - 1.0;
    return;
}
#endif

/*  Test problem test8
    # of real variables = 5
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 6
    */

#ifdef test8
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = 5.3578547 * x[2]*x[2] + 0.8356891*x[0]*x[4] + 37.293239*x[0]-40792.141;
    g[0] = 1.0 + (0.0056858*x[1]*x[4]+0.0006262*x[0]*x[3]-0.0022053*x[2]*x[4])/85.334407;
    g[1] = 1.0 - (0.0056858*x[1]*x[4]+0.0006262*x[0]*x[3]-0.0022053*x[2]*x[4])/6.665593;
    g[2] = -1.0 + (0.0071317*x[1]*x[4]+0.0029955*x[0]*x[1]+0.0021813*x[2]*x[2])/9.48751;
    g[3] = 1.0 - (0.0071317*x[1]*x[4]+0.0029955*x[0]*x[1]+0.0021813*x[2]*x[2])/29.48751;
    g[4] = -1.0 + (0.0047026*x[2]*x[4] + 0.0012547*x[0]*x[2]+0.0019085*x[2]*x[3])/10.699039;
    g[5] = 1.0 - (0.0047026*x[2]*x[4] + 0.0012547*x[0]*x[2]+0.0019085*x[2]*x[3])/15.699039;
     return;
}
#endif

/*  Test problem test9
    # of real variables = n
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 2
    */

#ifdef test9
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    int i;
    double sumcos4, sum2, prodcos2;
    double prod, sum;
    sumcos4 = sum2 = 0.0; prodcos2 = 1.0;
    for (i=0; i<nreal; i++)
    {
        sumcos4 +=  pow(cos(x[i]),4.0);
        prodcos2 *= pow(cos(x[i]),2.0);
        sum2 += (i+1.0)*x[i]*x[i];
    }
    obj[0] = fabs((sumcos4-2.0*prodcos2))/sqrt(sum2);
    prod = 1.0; sum = 0.0;
    for (i=0; i<nreal; i++)
    {
        prod *= x[i];
        sum  += x[i];
    }
    g[0] = prod/0.75 - 1.0;
    g[1] = 1.0 - sum/(7.5*nreal);
    return;
}
#endif

/*  Test problem test10
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 1
    */

#ifdef test10
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    double prod, sum;
    int i;
    prod= 1.0; sum=0.0;
    for (i=0; i<nreal; i++)
    {
        prod *= x[i];
        sum  += x[i]*x[i];
    }
    obj[0] = -1.0 * pow(nreal, nreal/2.0) * prod;
    if ((fabs(sum-1.0) <= EPS) && (obj[0] < -1.0))
    {
        obj[0] = -1.0;
    }
    sum = 0.0;
    for (i=0; i<nreal; i++)
    {
        sum  += x[i]*x[i];
    }
    g[0] = EPS - fabs(sum - 1.0);
    return;
}
#endif

/*  Test problem test11
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 3
    */

#ifdef test11
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = -1.0*(sin(2*PI*x[0]) * sin(2*PI*x[1]) / (square(x[0]) * x[0] * (x[0]+x[1])));
    g[0] = -1.0 * (x[0] * x[0] - x[1] + 1.0);
    g[1] = -1.0 * (1.0 - x[0] + (x[1]-4)*(x[1]-4.0));
    g[2] = -1.0 * (1.0 - (x[0]- (x[1]-1)*(x[1]-1.0))/0.1);
    return;
}
#endif

/*  Test problem test12
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 2
    */

#ifdef test12
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    obj[0] = (x[0]-10.0)*(x[0]-10.0)*(x[0]-10.0) + (x[1]-20.0)*(x[1]-20.0)*(x[1]-20.0);
    g[0] = (square(x[0]-5.0) + square(x[1]-5.0))/100.0 - 1.0;
    g[1] = 1.0 - (square(x[0]-6.0) + square(x[1]-5.0))/82.81;
    return;
}
#endif

/*  Test problem beam
    # of real variables = 4
    # of bin variables = 0
    # of objectives = 1
    # of constraints = 5
    */

#ifdef beam
void test_problem (double *x, double *b, int **gene, double *obj, double *g)
{
    double td, sd, force, emod, gmod, t1, t2, ctheta, ei, alpha, rad, mom, jmod, deltax;
    double tx, sx, pc;
    obj[0] = (1.0+0.37*0.283)*x[0]*x[0]*x[1] + 0.17*0.283 * x[2]*x[3]*(14.0+x[1]);
    td = 13600.0; sd = 30000.0; gmod = 12.0e06; emod = 30.0e06; force = 6000.0;
    mom = force * (14.0 + x[1]/2.0);
    rad = sqrt(x[1]*x[1]/4.0 + square(0.5*(x[2]+x[0])));
    jmod = 2.0 * (x[0]*x[1]*(x[1]*x[1]/12.0 + square(0.5*(x[0]+x[2])))/sqrt(2.0));
    t1 = force/(sqrt(2.0) * x[0] * x[1]);
    t2 = mom * rad / jmod;
    ctheta = x[1]/(2.0 * rad);
    tx = sqrt(square(t1) + square(t2) + 2.0*t1*t2*ctheta);
    sx = 6.0 * force * 14.0 / (x[3] * x[2] * x[2]);
    ei = emod * x[2] * x[3] * x[3] * x[3] / 12.0;
    alpha = gmod * x[2] * x[3] * x[3] * x[3] / 3.0;
    pc = 4.013 * sqrt(ei * alpha) * (1.0 - x[2] * sqrt(ei / alpha) / (2.0*14.0)) / (14.0 * 14.0);
    deltax = 4.0 * force * 14.0 * 14.0 * 14.0 / (emod * x[2]*x[2]*x[2] * x[3]);
    g[0] = 1.0 - tx/td;
    g[1] = 1.0 - sx/sd;
    g[2] = 1.0 - x[0]/x[3];
    g[3] = pc/force - 1.0;
    g[4] = 1.0 - deltax/0.25;
    return;
}
#endif
