
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "rand.h"
#include "qpoly.h"
#include "bvar.h"

qpoly_t *qpoly_new_zero(int n)
{
    qpoly_t *qpoly;
    qpoly = malloc(sizeof(qpoly_t));
    qpoly->n = n;
    qpoly->a = calloc(n + 1, sizeof(bvar_t)); /* 'n + 1' avoids to allocate 0 bytes when 'n = 0' */
    qpoly->b = 0;
    
    return qpoly;
}

qpoly_t *qpoly_new_random_linear(int n)
{
    qpoly_t *qpoly = qpoly_new_zero(n);

    int i;
    for (i = 0; i < n; i++)
    {
        if (rand_bool())
        {
            qpoly->a[i] |= (1 << i);
        }
    }

    qpoly->b = rand_bool();

    return qpoly;
}

qpoly_t *qpoly_new_random(int n)
{
    qpoly_t *qpoly = qpoly_new_zero(n);
    
    int i, j;        
    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {
            if (rand_bool())
            {
                qpoly->a[i] |= (1 << j);
            }
        }
    }

    qpoly->b = rand_bool();
    
    return qpoly;
}

void qpoly_free(qpoly_t *qpoly)
{
    free(qpoly->a);
    free(qpoly);
}

qpoly_t *qpoly_copy(qpoly_t *qpoly)
{
    qpoly_t *copy;
    
    copy = qpoly_new_zero(qpoly->n);
    memcpy(copy->a, qpoly->a, qpoly->n * sizeof(bvar_t));
    copy->b = qpoly->b;
    
    return copy;
}

void qpoly_print(qpoly_t *qpoly)
{
    int i, j, p;
    
    for (i = 0, p = 0; i < qpoly->n; i++)
    {
        for (j = i; j < qpoly->n; j++)
        {
            if (qpoly->a[i] & (1 << j))
            {
                if (p)
                {
                    printf(" + ");
                }
                
                if (i == j)
                {
                    printf("x_%d", i + 1);
                } 
                else
                {
                    printf("x_%d x_%d", i + 1, j + 1);
                }
                p = 1;
            }
        }
    }
    
    if (!p)
    {
        printf("%d", qpoly->b);
    }
    else if(qpoly->b)
    {
        printf(" + 1");
    }
    
    printf("\n");
}

bool qpoly_eval(qpoly_t *qpoly, bvar_t x)
{
    bvar_t y;
    int r, i;
    for (i = 0, r = qpoly->b, y = x; i < qpoly->n; i++, y >>= 1)
    {
        if (y & 1)
        {
            r ^= __builtin_parityl(qpoly->a[i] & x);
        }
    }
    
    return r;
}

void qpoly_add(qpoly_t *qpoly1, qpoly_t *qpoly2)
{
    int i;
    for (i = 0; i < qpoly1->n; i++)
    {
        qpoly1->a[i] ^= qpoly2->a[i];
    }
    qpoly1->b ^= qpoly2->b;
}

void qpoly_subs(qpoly_t *qpoly, bvar_t z, int k)
{
    int n = qpoly->n;
    
    int i;
    for (i = 0; i < n; i++)
    {
        if (__builtin_parity((qpoly->a[i] >> (n - k)) & z))
        {
            if (i < n - k)
            {
                qpoly->a[i] ^= (1 << i);
            }
            else
            {
                if (z & (1 << (i - n + k)))
                {
                    qpoly->b ^= 1;
                }
            }
        }
        
        qpoly->a[i] &= (1 << (n - k)) - 1;
    }
    
    
    qpoly->n = n - k;
}

bfunc_t *qpoly_to_bfunc(qpoly_t *qpoly)
{
    bfunc_t *bfunc;
    
    bfunc = bfunc_new(qpoly->n);

    int i, j;
    for (i = 0; i < qpoly->n; i++)
    {
        for (j = i; j < qpoly->n; j++)
        {
            if ((qpoly->a[i] >> j) & 1)
            {
                bvar_t mon;
                
                mon = ((bvar_t)1 << i) | ((bvar_t)1 << j);
                
                bfunc_set(bfunc, mon, 1);
            }
        }
    }
    
    if (qpoly->b)
    {
        bfunc_set(bfunc, 0, 1);
    }
    
    return bfunc;
}

bool qpoly_costant_term(qpoly_t *qpoly)
{
    return qpoly->b;
}

