#ifndef BINOMIALS_H
#define BINOMIALS_H

#include <stdint.h>

/* Return the binomial coefficient (n, k).
 * (Assuming 0 <= n, k <= 60.) */
uint64_t binomial(int n, int k);

/* Return sum_{i = 0}^k binomial(n, k).
 * (Assuming 0 <= n, k <= 60.) */
uint64_t sum_binomials(int n, int k);

#endif
