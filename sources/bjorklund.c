
#include <math.h>
#include <stdlib.h>
#include "bfunc.h"
#include "bjorklund.h"

bool bjorklund_parity(qsyst_t *qsyst, double lambda, int c)
{
    int n, n1, n2, d, s;
    
    n = qsyst_num_vars(qsyst);
    n1 = (int)floor(lambda * n);
    n2 = n - n1;

    if (n1 == 0)
    {
        return qsyst_count(qsyst) % 2;
    }

    d = n1 + 4;
    s = 48 * n + 1;
    
    if (c) 
    {
        s = c;
    }
    
    bfunc_t **bfunc;
    
    bfunc = malloc(s * sizeof(bfunc_t*));
    
    int k;
    for (k = 0; k < s; k++)
    {
        bfunc[k] = bfunc_new(n2);
        
        qsyst_t *R = qsyst_rand_lin_comb(qsyst, n1 + 2);
        
        bvar_t z = 0; int h = 0;
        do
        {
            qsyst_t *Q = qsyst_copy(R);
            qsyst_subs(Q, z, n2);
            
            bfunc_set(bfunc[k], z, bjorklund_parity(Q, lambda, c));
            
            qsyst_free(Q);
        }
        while(next_subset(&z, &h, n2, d));
        
        bfunc_restricted_zeta_transform(bfunc[k], d);
        bfunc_zeta_transform(bfunc[k]);

        qsyst_free(R);
    }
    
    int I = 0;
    
    bvar_t z;
    for (z = 0; z < ((bvar_t)1 << n2); z++)
    {
        int c = 0;
        for (k = 0; k < s; k++)
        {
            c += bfunc_get(bfunc[k], z);
        }

        if (2 * c > s)
        {
            I ^= 1;
        }
    }
    
    for (k = 0; k < s; k++)
    {
        bfunc_free(bfunc[k]);
    }
    
    free(bfunc);
    
    return I;

}
