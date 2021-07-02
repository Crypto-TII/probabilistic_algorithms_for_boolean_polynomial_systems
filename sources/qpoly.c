
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "rand.h"
#include "qpoly.h"
#include "bpoly.h"
#include "bvar.h"


/* Structure for a quadratic (Boolean) polynomial
 * 
 * Q(x) = sum_{1 <= i <= j < n} a_{i,j} x_i x_j + b,  a_{i,j}, b in F_2
 * 
 * in the variables x_1, ..., x_n (note that x_i^2 = x_i).
 */
struct qpoly_s
{
    int     n; /* Number of variables. */
    bvar_t *a; /* (a[i] << j) & 1 is equal to a_{i,j}, if i <= j, 
                                       and to       0, otherwise. */
    bool    b; /* Constant term. */
};

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

bpoly_t *qpoly_to_bpoly(qpoly_t *qpoly)
{
    bpoly_t *bpoly = bpoly_new_zero();
    
    if (qpoly->b)
    {
        bpoly_add_mon(bpoly, 0);
    }
    
    int i, j;
    for (i = 0; i < qpoly->n; i++)
    {
        for (j = i; j < qpoly->n; j++)
        {
            if ((qpoly->a[i] >> j) & 1)
            {
                bvar_t x = ((bvar_t)1 << i) | ((bvar_t)1 << j);
                
                bpoly_add_mon(bpoly, x);
            }
        }
    }

    return bpoly;
}

bfunc_t *qpoly_to_bfunc(qpoly_t *qpoly)
{
    bpoly_t *bpoly = qpoly_to_bpoly(qpoly);
    bfunc_t *bfunc = bpoly_monomials(bpoly, qpoly->n);
    bpoly_free(bpoly);
    return bfunc;
}

bool qpoly_costant_term(qpoly_t *qpoly)
{
    return qpoly->b;
}

