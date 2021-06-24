/* A custom doubly linked list implemenation */

#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "rand.h"

/* Insert an element X into the list at location specified by NODE */
void insert (list *node, int x)
{
    list *temp;
    if (node==NULL)
    {
        error("Asked to enter after a NULL pointer.");
    }
    temp = (list *)malloc(sizeof(list));
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}

/* Delte the node NODE from the list */
list* del (list *node)
{
    list *temp;
    if (node==NULL)
    {
        error("Asked to delete a NULL pointer.");
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child!=NULL)
    {
        temp->child->parent = temp;
    }
    free (node);
    return (temp);
}
