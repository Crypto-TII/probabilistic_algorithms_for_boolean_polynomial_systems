
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "rbfunc.h"

typedef uint64_t rbchunk;

const uint64_t rbchunk_size = 8 * sizeof(rbchunk);

struct rbfunc_s
{
    int        n; /* Number of variables. */
    int        w; /* upper bound weigth. */
    bvar_t  nvalues; /* Number of values. */
    rbchunk *data; /* Truth values or monomial coefficients. */
};

rbfunc_t *rbfunc_new(int n, int w)
{
    rbfunc_t *rbfunc;

    rbfunc = malloc(sizeof(rbfunc_t));
    rbfunc->n = n;
    rbfunc->w = (w < n) ? w : n;
    rbfunc->nvalues = bvar_number_of_values(n, rbfunc->w);
    rbfunc->data = calloc(rbfunc->nvalues / rbchunk_size + 1, sizeof(rbchunk));

    return rbfunc;
}

/* Replace 'rbfunc1' with 'rbfunc1 + bfunc2 * mon' */
void rbfunc_mul_add(rbfunc_t *rbfunc1, rbfunc_t *rbfunc2, bvar_t mon)
{
    bvar_t x = 0; int h = 0;
    do
    {
        rbfunc_xor(rbfunc1, x | mon, rbfunc_get(rbfunc2, x));
    }
    while (next_subset(&x, &h, rbfunc2->n, rbfunc2->w));
}

/* Replace 'rbfunc1' with 'rbfunc2 * (qpoly + 1)' */
void rbfunc_mul_qpoly(rbfunc_t *rbfunc1, rbfunc_t *rbfunc2, qpoly_t *qpoly)
{
    int i, j;
    for (i = 0; i < qpoly->n; i++)
    {
        for (j = i; j < qpoly->n; j++)
        {
            if (qpoly->a[i] & (1 << j))
            {
                bvar_t mon;
                
                mon = ((bvar_t)1 << i) | ((bvar_t)1 << j);

                rbfunc_mul_add(rbfunc1, rbfunc2, mon);
            }
        }
    }
    
    if (!qpoly->b)
    {
        rbfunc_mul_add(rbfunc1, rbfunc2, 0);
    }
}

rbfunc_t *rbfunc_new_characteristic(qsyst_t *qsyst)
{
    int n, m, w;
    rbfunc_t *rbfunc2;
    
    n = qsyst_num_vars(qsyst);
    m = qsyst_num_eqs(qsyst);
    w = (2*m < n) ? (2*m) : n;

    rbfunc2 = rbfunc_new(n, w);
    rbfunc_set(rbfunc2, 0, 1);
    
    int i;
    for (i = 0; i < m; i++)
    {
        rbfunc_t *rbfunc1, *rbfunc3;
        
        rbfunc1 = rbfunc_new(n, w);
        
        rbfunc_mul_qpoly(rbfunc1, rbfunc2, qsyst_equ(qsyst, i));
        
        rbfunc3 = rbfunc2;
        rbfunc2 = rbfunc1;
        rbfunc_free(rbfunc3);
    }
    
    return rbfunc2;
}

void rbfunc_free(rbfunc_t *rbfunc)
{
    free(rbfunc->data);
    free(rbfunc);
}

bool rbfunc_get(rbfunc_t *rbfunc, bvar_t x)
{
    bvar_t x_index = bvar_get_index(x, rbfunc->n);
    bvar_t number_of_values = rbfunc->nvalues;
    
    if (number_of_values < x_index)
    {
        return 0;
    }
    else
    {
       return (rbfunc)->data[(x_index) / rbchunk_size] & ((rbchunk)1 << ((x_index) % rbchunk_size));
    }

}

void rbfunc_set(rbfunc_t *rbfunc, bvar_t x, bool y)
{
    bvar_t x_index = bvar_get_index(x, rbfunc->n);
    bvar_t number_of_values = rbfunc->nvalues;
    if (x_index < number_of_values)
    {
        if (y)
        {
            (rbfunc)->data[(x_index) / rbchunk_size] |= ((rbchunk)1 << ((x_index) % rbchunk_size));
        }
        else
        {
            (rbfunc)->data[(x_index) / rbchunk_size] &= ~((rbchunk)1 << ((x_index) % rbchunk_size));
        }
    }
}

rbfunc_t *rbfunc_copy(rbfunc_t *rbfunc)
{
    rbfunc_t *copy = rbfunc_new(rbfunc->n, rbfunc->w);
    
    bvar_t x = 0; int h = 0;
    do
    {
        rbfunc_set(copy, x, rbfunc_get(rbfunc, x));
    }
    while (next_subset(&x, &h, rbfunc->n, rbfunc ->w));
    
    return copy;
}

void rbfunc_add(rbfunc_t *rbfunc1, rbfunc_t *rbfunc2)
{
    bvar_t i;
    for (i = 0; i < rbfunc1->nvalues / rbchunk_size + 1; i++)
    {
        rbfunc1->data[i] ^= rbfunc2->data[i];
    }
}

void rbfunc_print(rbfunc_t *rbfunc)
{
    bvar_t x = 0; int h = 0;
    do
    {
        bvar_print_map(x, rbfunc->n, rbfunc_get(rbfunc, x));
    }
    while (next_subset(&x, &h, rbfunc->n, rbfunc ->w));
}


void rbfunc_xor(rbfunc_t *rbfunc, bvar_t x, bool y)
{
    bvar_t x_index = bvar_get_index(x, rbfunc->n);
    (rbfunc)->data[(x_index) / rbchunk_size] ^= ((y) ? ((rbchunk)1 << ((x_index) % rbchunk_size)) : 0);
}

void rbfunc_zeta_transform(rbfunc_t *rbfunc)
{
    int i;
    for (i = 0; i < rbfunc->n; i++)
    {
        bvar_t x = 0; int h = 0;
        do
        {
            if (x & ((bvar_t)1 << i))
            {
                rbfunc_xor(rbfunc, x, rbfunc_get(rbfunc, x ^ ((bvar_t)1 << i)));
            }
        }
        while (next_subset(&x, &h, rbfunc->n, rbfunc ->w));
    }
}

bfunc_t *rbfunc_to_bfunc(rbfunc_t *rbfunc)
{
    int n = rbfunc->n;
    int w = rbfunc->w;
    
    bfunc_t *bfunc = bfunc_new(n);
    
    bvar_t x = 0; int h = 0;
    do
    {
        bfunc_set(bfunc, x, rbfunc_get(rbfunc, x));
    }
    while (next_subset(&x, &h, n, w));
    
    return bfunc;
}
