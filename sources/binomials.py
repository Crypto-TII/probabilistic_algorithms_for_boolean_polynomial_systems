
# Run this script to generate 'binomials.c'

from sympy import binomial

f = open("binomials.c", "w")

f.write("""
#include "binomials.h"

uint64_t pre_binomial[] = {
""")

for n in range(0, 60 + 1):
    for k in range(0, n + 1):
        f.write("(uint64_t)" + str(binomial(n, k)) + ", ")
    f.write("\n")
f.write("};")

f.write("""

uint64_t pre_sum_binomial[] = {
""")

for n in range(0, 60 + 1):
    S = 0
    for k in range(0, n + 1):
        S += binomial(n, k)
        f.write("(uint64_t)" + str(S) + ", ")
    f.write("\n")
f.write("};")

f.write("""

uint64_t binomial(int n, int k)
{
    if (k > n)
    {
        return 0;
    }
    
    return pre_binomial[(n * (n + 1) / 2) + k];
}

uint64_t sum_binomials(int n, int k)
{
    if (k >= n)
    {
        return ((uint64_t)1 << n);
    }
    
    return pre_sum_binomial[(n * (n + 1) / 2) + k];
}

""")

f.close()
