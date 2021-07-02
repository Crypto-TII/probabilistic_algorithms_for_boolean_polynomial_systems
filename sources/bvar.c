
#include <stdio.h>
#include "binomials.h"
#include "bvar.h"

void bvar_print(bvar_t x, int n)
{
    char s[BVAR_MAX + 1];
    int i;
    for (i = n - 1; i >= 0; i--, x >>= 1)
    {
        s[i] = (x & 1) ? '1' : '0';
    }
    
    s[n] = '\0';
    
    printf("%s\n", s);    
}

void bvar_print_map(bvar_t x, int n, bool y)
{
    char s[BVAR_MAX + 1];
    int i;
    for (i = n - 1; i >= 0; i--, x >>= 1)
    {
        s[i] = (x & 1) ? '1' : '0';
    }
    
    s[n] = '\0';
    
    printf("%s -> %d\n", s, y);        
}

bool next_subset(bvar_t *x, int *h, int n, int k)
{
    if (*x == 0)
    {
        if (k > 0)
        {
            *x = 1;
            *h = 1;
            return 1;
        }
        else
        {
            return 0;
        }
    }

    bvar_t t = *x;
    
    /* Compute next subset as explained in Item 175 of HAKMEM * * * * */
    bvar_t smallest, ripple, ones;
    smallest = (*x) & (-*x);
    ripple = (*x) + smallest;
    ones = (*x) ^ ripple;
    ones = (ones >> 2) / smallest;
    
    *x = (ripple | ones) & (((bvar_t)1 << n) - 1);
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */
    
    if (*x < t)
    {
        (*h)++;
        
        if (*h <= k)
        {
            *x = ((bvar_t)1 << *h) - 1;
            return 1;
        }
        else
        {
            return 0;
        }
    }
    
    return 1;
}

bvar_t bvar_number_of_values(int n, int w)
{
    return sum_binomials(n, w);
}

bvar_t bvar_get_index(bvar_t b, int n)
{
    int b_weight = bvar_weight(b);
    bvar_t x = 0;
    int counter = 0;
    if (b)
    {
        int i;
        for (i = 0; i < n; ++i) 
        {
            if ((b >> i) & 1)
            {
                counter += 1;
                x += binomial(i, counter);
            }
        }
        
        return x + bvar_number_of_values(n, b_weight - 1);
    }
    else
    {
        return 0;
    }
}

int bvar_weight(bvar_t x)
{
    return __builtin_popcount(x);
}

bvar_t bvar_first_bits(bvar_t x, int n1){
    return x & ((1 << n1) - 1);
}

/* Returns last n - n1 bits of x, where x is an bvar object of lenght n */
bvar_t bvar_last_bits(bvar_t x, int n, int n1){
    return (x & ((1 << (n - n1)) - 1) << n1) >> n1;
}
