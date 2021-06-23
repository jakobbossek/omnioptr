#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"

void selection (population *old_pop, population *new_pop)
{
    int *a1, *a2;
    int temp;
    int i;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *)malloc(popsize*sizeof(int));
    a2 = (int *)malloc(popsize*sizeof(int));
    for (i=0; i<popsize; i++)
    {
        a1[i] = a2[i] = i;
    }
    for (i=0; i<popsize; i++)
    {
        rand = rnd (i, popsize-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rnd (i, popsize-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    for (i=0; i<popsize; i+=4)
    {
        parent1 = tournament (&old_pop->ind[a1[i]], &old_pop->ind[a1[i+1]]);
        parent2 = tournament (&old_pop->ind[a1[i+2]], &old_pop->ind[a1[i+3]]);
        crossover (parent1, parent2, &new_pop->ind[i], &new_pop->ind[i+1]);
        parent1 = tournament (&old_pop->ind[a2[i]], &old_pop->ind[a2[i+1]]);
        parent2 = tournament (&old_pop->ind[a2[i+2]], &old_pop->ind[a2[i+3]]);
        crossover (parent1, parent2, &new_pop->ind[i+2], &new_pop->ind[i+3]);
    }
    free (a1);
    free (a2);
    return;
}

void restricted_selection (population *old_pop, population *new_pop)
{
    int *a1, *a2;
    list *n1, *n2;
    list *temp1, *temp2;
    list *temp11;
    double *min_realvar, *max_realvar;
    double *min_binvar, *max_binvar;
    int temp;
    int i, j;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *)malloc(popsize*sizeof(int));
    a2 = (int *)malloc(popsize*sizeof(int));
    n1 = (list *)malloc(sizeof(list));
    n1->index = -1;
    n1->parent = NULL;
    n1->child = NULL;
    n2 = (list *)malloc(sizeof(list));
    n2->index = -1;
    n2->parent = NULL;
    n2->child = NULL;
    min_realvar = NULL;
    max_realvar = NULL;
    min_binvar = NULL;
    max_binvar = NULL;
    if (nreal!=0)
    {
        min_realvar = (double *)malloc(nreal*sizeof(double));
        max_realvar = (double *)malloc(nreal*sizeof(double));
        for (i=0; i<nreal; i++)
        {
            min_realvar[i] = INF;
            max_realvar[i] = -INF;
        }
        for (i=0; i<nreal; i++)
        {
            for (j=0; j<popsize; j++)
            {
                if (min_realvar[i] > old_pop->ind[j].xreal[i])
                {
                    min_realvar[i] = old_pop->ind[j].xreal[i];
                }
                if (max_realvar[i] < old_pop->ind[j].xreal[i])
                {
                    max_realvar[i] = old_pop->ind[j].xreal[i];
                }
            }
        }
    }
    if (nbin!=0)
    {
        min_binvar = (double *)malloc(nbin*sizeof(double));
        max_binvar = (double *)malloc(nbin*sizeof(double));
        for (i=0; i<nbin; i++)
        {
            min_binvar[i] = INF;
            max_binvar[i] = -INF;
        }
        for (i=0; i<nbin; i++)
        {
            for (j=0; j<popsize; j++)
            {
                if (min_binvar[i] > old_pop->ind[j].xbin[i])
                {
                    min_binvar[i] = old_pop->ind[j].xbin[i];
                }
                if (max_binvar[i] < old_pop->ind[j].xbin[i])
                {
                    max_binvar[i] = old_pop->ind[j].xbin[i];
                }
            }
        }
    }
    for (i=0; i<popsize; i++)
    {
        a1[i] = a2[i] = i;
    }
    for (i=0; i<popsize; i++)
    {
        rand = rnd (i, popsize-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rnd (i, popsize-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    temp1 = n1;
    temp2 = n2;
    for (i=0; i<popsize; i++)
    {
        insert(temp1,a1[i]);
        insert(temp2,a2[i]);
        temp1 = temp1->child;
        temp2 = temp2->child;
    }
    i=0;
    do
    {
        temp1 = n1->child;
        temp11 = search_nearest (old_pop, temp1, n1, min_realvar, max_realvar, min_binvar, max_binvar);
        parent1 = tournament (&old_pop->ind[temp1->index], &old_pop->ind[temp11->index]);
        temp1 = del(temp1);
        temp11 = del(temp11);
        temp1 = n1->child;
        temp11 = search_nearest (old_pop, temp1, n1, min_realvar, max_realvar, min_binvar, max_binvar);
        parent2 = tournament (&old_pop->ind[temp1->index], &old_pop->ind[temp11->index]);
        temp1 = del(temp1);
        temp11 = del(temp11);
        crossover (parent1, parent2, &new_pop->ind[i], &new_pop->ind[i+1]);

        temp1 = n2->child;
        temp11 = search_nearest (old_pop, temp1, n2, min_realvar, max_realvar, min_binvar, max_binvar);
        parent1 = tournament (&old_pop->ind[temp1->index], &old_pop->ind[temp11->index]);
        temp1 = del(temp1);
        temp11 = del(temp11);
        temp1 = n2->child;
        temp11 = search_nearest (old_pop, temp1, n2, min_realvar, max_realvar, min_binvar, max_binvar);
        parent2 = tournament (&old_pop->ind[temp1->index], &old_pop->ind[temp11->index]);
        temp1 = del(temp1);
        temp11 = del(temp11);
        crossover (parent1, parent2, &new_pop->ind[i+2], &new_pop->ind[i+3]);
        i+=4;
    }
    while (temp1->child != NULL);
	if (nreal!=0)
	{
		free (min_realvar);
		free (max_realvar);
	}
	if (nbin!=0)
	{
		free (min_binvar);
		free (max_binvar);
	}
    free (a1);
    free (a2);
	free (n1);
	free (n2);
    return;
}

individual* tournament (individual *ind1, individual *ind2)
{
    int flag;
    flag = check_dominance (ind1, ind2);
    if (flag==1)
    {
        return (ind1);
    }
    if (flag==-1)
    {
        return (ind2);
    }
    if (ind1->crowd_dist > ind2->crowd_dist)
    {
        return(ind1);
    }
    if (ind2->crowd_dist > ind1->crowd_dist)
    {
        return(ind2);
    }
    if ((randomperc()) <= 0.5)
    {
        return(ind1);
    }
    else
    {
        return(ind2);
    }
}

list* search_nearest (population *pop, list *temp1, list *n1, double *min_realvar, double *max_realvar, double *min_binvar, double *max_binvar)
{
    list *temp11;
    list *temp;
    double min_dist;
    double dist;
    if (n1->child->child==NULL)
    {
        printf("\n Error in restricted mating, hence exiting \n");
        exit(1);
    }
    if (n1->child->child->child==NULL)
    {
        return (n1->child->child);
    }
    min_dist = INF;
    temp11=n1->child->child;
    temp = temp11;
    do
    {
        dist = calc_distance (pop, temp1, temp, min_realvar, max_realvar, min_binvar, max_binvar);
        if (dist < min_dist)
        {
            temp11 = temp;
            min_dist = dist;
        }
        temp = temp->child;
    }
    while (temp!=NULL);
    return (temp11);
}

double calc_distance (population *pop, list *temp1, list *temp, double *min_realvar, double *max_realvar, double *min_binvar, double *max_binvar)
{
    double distance;
    double res;
    int i;
    distance = 0.0;
    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            if (min_realvar[i] == max_realvar[i])
            {
                distance += 0.0;
            }
            else
            {
                res = (pop->ind[temp1->index].xreal[i]-pop->ind[temp->index].xreal[i])/(max_realvar[i]-min_realvar[i]);
                distance += pow(res,2.0);
            }
        }
    }
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            if (min_binvar[i] == max_binvar[i])
            {
                distance += 0.0;
            }
            else
            {
                res = (pop->ind[temp1->index].xbin[i]-pop->ind[temp->index].xbin[i])/(max_binvar[i]-min_binvar[i]);
                distance += pow(res,2.0);
            }
        }
    }
    distance = sqrt (distance);
    return (distance);
}
