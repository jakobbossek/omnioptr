#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"

/* Routine for modified epsilon dominance criterion
 * It will return the following values
 * 1 if a dominates b
 * -1 if b dominates a
 * 0 if both a and b are non-dominated
 */
int check_loose_dominance (individual *a, individual *b)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;
    if (a->constr_violation<0 && b->constr_violation<0)
    {
        if (a->constr_violation > b->constr_violation)
        {
            return (1);
        }
        else
        {
            if (a->constr_violation < b->constr_violation)
            {
                return (-1);
            }
            else
            {
                return (0);
            }
        }
    }
    else
    {
        if (a->constr_violation < 0 && b->constr_violation == 0)
        {
            return (-1);
        }
        else
        {
            if (a->constr_violation == 0 && b->constr_violation <0)
            {
                return (1);
            }
            else
            {
                for (i=0; i<nobj; i++)
                {
                    if (a->obj[i] < b->obj[i])
                    {
                        flag1 = 1;

                    }
                    else
                    {
                        if (a->obj[i] > b->obj[i])
                        {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1==1 && flag2==0)
                {
                    flag1 = 0;
                    flag2 = 0;
                    for (i=0; i<nobj; i++)
                    {
                        if (a->obj[i]+(delta*epsilon[i]) < b->obj[i])
                        {
                            flag1 = 1;
                        }
                    }
                    if (flag1 == 1)
                    {
                        return (1);
                    }
                    else
                    {
                        return (0);
                    }
                }
                else
                {
                    if (flag1==0 && flag2==1)
                    {
                        flag1 = 0;
                        flag2 = 0;
                        for (i=0; i<nobj; i++)
                        {
                            if (a->obj[i] > b->obj[i]+(delta*epsilon[i]))
                            {
                                flag2 = 1;
                            }
                        }
                        if (flag2 == 1)
                        {
                            return (-1);
                        }
                        else
                        {
                            return (0);
                        }
                    }
                    else
                    {
                        return (0);
                    }
                }
            }
        }
    }
}

/* Routine for usual non-domination checking
 * It will return the following values
 * 1 if a dominates b
 * -1 if b dominates a
 * 0 if both a and b are non-dominated
*/
int check_dominance (individual *a, individual *b)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;
    if (a->constr_violation<0 && b->constr_violation<0)
    {
        if (a->constr_violation > b->constr_violation)
        {
            return (1);
        }
        else
        {
            if (a->constr_violation < b->constr_violation)
            {
                return (-1);
            }
            else
            {
                return (0);
            }
        }
    }
    else
    {
        if (a->constr_violation < 0 && b->constr_violation == 0)
        {
            return (-1);
        }
        else
        {
            if (a->constr_violation == 0 && b->constr_violation <0)
            {
                return (1);
            }
            else
            {
                for (i=0; i<nobj; i++)
                {
                    if (a->obj[i] < b->obj[i])
                    {
                        flag1 = 1;

                    }
                    else
                    {
                        if (a->obj[i] > b->obj[i])
                        {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1==1 && flag2==0)
                {
                    return (1);
                }
                else
                {
                    if (flag1==0 && flag2==1)
                    {
                        return (-1);
                    }
                    else
                    {
                        return (0);
                    }
                }
            }
        }
    }
}

void define_epsilon (population *pop, int size, double *eps)
{
    int i, j;
    for (j=0; j<nobj; j++)
    {
        min_obj[j] = pop->ind[0].obj[j];
        max_obj[j] = pop->ind[0].obj[j];
    }
    for (j=0; j<nobj; j++)
    {
        for (i=0; i<size; i++)
        {
            if (pop->ind[i].obj[j] < min_obj[j])
            {
                min_obj[j] = pop->ind[i].obj[j];
            }
            if (pop->ind[i].obj[j] > max_obj[j])
            {
                max_obj[j] = pop->ind[i].obj[j];
            }
        }
    }
    for (j=0; j<nobj; j++)
    {
        eps[j] = max_obj[j] - min_obj[j];
    }
    return;
}
