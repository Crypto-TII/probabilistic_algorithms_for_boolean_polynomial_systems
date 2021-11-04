
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "rand.h"
#include "bfunc.h"
#include "rbfunc.h"
#include "lokshtanov.h"

bool lokshtanov_consistency(qsyst_t *qsyst, double delta, int number_of_iterations){
    
    int n, n1, s;
    
    n = qsyst_num_vars(qsyst);
    n1 = (int)floor(delta * n);
    s = 100 * n;
    
    if (number_of_iterations)
    {
        s = number_of_iterations;
    }
    
    if(n == 0)
    {
        return qsyst_are_all_zero(qsyst);
    }
    
    int *score = calloc(1 << (n - n1), sizeof(int));
    
    int w = (2*(n1 + 2) < n - n1) ? (2*(n1 + 2)) : (n - n1);
    
    int i;
    for (i = 0; i < s; i++) 
    {    
        rbfunc_t *R = rbfunc_new(n - n1, w);
        
        bvar_t z;
        for (z = 0; z < ((bvar_t)1 << n1); z++) 
        {
            if (rand_bool()) 
            {
                qsyst_t *subs, *rlc;
                rbfunc_t *F;
                
                subs = qsyst_copy(qsyst);
                qsyst_subs(subs, z, n1);

                rlc = qsyst_rand_lin_comb(subs, n1 + 2);
                F = rbfunc_new_characteristic(rlc);

                rbfunc_add(R, F);
                
                qsyst_free(subs);
                qsyst_free(rlc);
                rbfunc_free(F);
            }
        }
        
        bfunc_t *bfunc;
        
        bfunc = rbfunc_to_bfunc(R);
        rbfunc_free(R);
        bfunc_zeta_transform(bfunc);
        
        bvar_t x;
        for (x = 0; x < ((bvar_t)1 << (n - n1)); x++) 
        {
            if (bfunc_get(bfunc, x))
            {
                score[x] += 1;
            }
        }
        
        bfunc_free(bfunc);
    }
    
    int threshold = 0.4 * s;
    bool cons = 0;
    
    bvar_t x;
    for (x = 0; x < ((bvar_t)1 << (n - n1)); x++)
    {
        if (score[x] > threshold)
        {
            cons = 1;
            break;
        }
    }
    
    free(score);
    
    return cons;
}
