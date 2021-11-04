
#include <stdlib.h>
#include <math.h>
#include "dinur.h"
#include "bfunc.h"

bool dinur_parity(qsyst_t *qsyst, double kappa, double lambda, int number_of_iterations)
{
    int n, n1;
    
    n  = qsyst_num_vars(qsyst);
    n1 = (int)floor(kappa * n);


    bfunc_t *V = dinur_multiparity(qsyst, n1, n - n1, lambda, number_of_iterations);
    
    bool p = 0;
    
    bvar_t y;
    for (y = 0; y < ((bvar_t)1 << (n - n1)); y++)
    {
        p ^= bfunc_get(V, y);
    }
    
    bfunc_free(V);
    
    return p;
}

bfunc_t *dinur_multiparity(qsyst_t *qsyst, int n1, int w, double lambda, int number_of_iterations)
{
    int n, n2, t;

    n = qsyst_num_vars(qsyst);
    n2 = (int)floor(n1 - lambda * n);
    t = 48 * n + 1;
    
    if (number_of_iterations)
    {
        t = number_of_iterations;
    }
    
    if (n2 <= 0)
    {
        return dinur_bruteforce_multiparity(qsyst, n1, w);
    }
    
    bfunc_t **score = malloc(t * sizeof(bfunc_t*));
    
    int k;
    for (k = 0; k < t; k++)
    {
        qsyst_t *rlc = qsyst_rand_lin_comb(qsyst, n2 + 2);
        
        bfunc_t *Vk = dinur_multiparity(rlc, n2, n2 + 4, lambda, number_of_iterations);
        
        bfunc_restricted_zeta_transform(Vk, n2 + 4);
        bfunc_partially_restricted_zeta_transform(Vk, n1 - n2, w);
        
        score[k] = Vk;
        
        qsyst_free(rlc);
    }
    
    bfunc_t *V = bfunc_new(n - n1);
    
    bvar_t y = 0; int h = 0;
    do
    {
        bvar_t u;
        for (u = 0; u < ((bvar_t)1 << (n1 - n2)); u++)
        {
            int c = 0;
            
            int k;
            for (k = 0; k < t; k++)
            {
                if (bfunc_get(score[k], y | (u << (n - n1))))
                {
                    c++;
                }
            }
            
            if (2 * c > t)
            {
                bfunc_add(V, y, 1);
            }
        }
    }
    while (next_subset(&y, &h, n - n1, w));

    for (k = 0; k < t; k++)
    {
        bfunc_free(score[k]);
    }

    free(score);
    
    return V;
}

bfunc_t *dinur_bruteforce_multiparity(qsyst_t * qsyst, int n1, int w)
{
    int m, n;
    
    m = qsyst_num_eqs(qsyst);
    n = qsyst_num_vars(qsyst);
    
    bfunc_t *evals = bfunc_new(n);
    
    bvar_t x;
    for (x = 0; x < ((bvar_t)1 << n); x++)
    {
        bfunc_set(evals, x, 1);
    }
    
    int j;
    for (j = 0; j < m; j++)
    {
        bfunc_t *bfunc;
        
        bfunc = qpoly_to_bfunc(qsyst_equ(qsyst, j));
        bfunc_add(bfunc, 0, 1);
        
        bfunc_partially_restricted_zeta_transform(bfunc, n1, w);
        
        bfunc_and(evals, bfunc);
        
        bfunc_free(bfunc);
    }
    
    bfunc_t *V = bfunc_new(n - n1);

    bvar_t y = 0; int h = 0;
    do
    {
        bvar_t z;
        for (z = 0; z < ((bvar_t)1 << n1); z++)
        {
            bfunc_add(V, y, bfunc_get(evals, y | (z << (n - n1))));
        }
    }
    while (next_subset(&y, &h, n - n1, w));
    
    bfunc_free(evals);

    return V;
}
