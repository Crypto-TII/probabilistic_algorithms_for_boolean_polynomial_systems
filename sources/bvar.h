#ifndef BVAR_H
#define BVAR_H

#include <stdbool.h>
#include <stdint.h>

/* Type for a vector of Boolean variables.
 * For example, 
 * 'bvar_t x = b1011' means 'x_1 = 1, x_2 = 1, x_3 = 0, x_4 = 1. */
typedef uint64_t bvar_t;

/* Maximum number of Boolean variables that can be stored
 * (and work with!) in bvar_t. */
const static int BVAR_MAX = 8 * sizeof(bvar_t) - 1;

/* Print the 'n' Boolean variables 'x'.
 * (For testing/debugging.) */
void bvar_print(bvar_t x, int n);

/* Print the 'n' Boolean variables 'x' and their associated value 'y'.
 * (For testing/debugging.) */
void bvar_print_map(bvar_t x, int n, bool y);

/* Execute 'stuff(x, h)' over all 'n'-bits binary words 
 * of Hamming weight at most 'k'
 * 
 * bvar_t x = 0; int h = 0;
 * do
 * {
 *     stuff(x, h);
 * }
 * while (next_subset(&x, &h, n, k));
 * 
 * NOTE: It must be 'n < BVAR_MAX'.
 */ 
bool next_subset(bvar_t *x, int *h, int n, int k);

/* Returns the number of values in a restricted boolean function object of 'n' variables and bound weight w. */
bvar_t bvar_number_of_values(int n, int w);

/* Return sum_{i=0}^w binomial(n, i), that is, the number 
 * of binary words of length 'n' and Hamming weigth <= 'w'. */
bvar_t bvar_sum_binomials(int n, int w);

/*output the the index of b in the set of integers of lenght n and weight <= w*/
bvar_t bvar_get_index(bvar_t b, int n);

/* Returns the Hamming weight of the Boolean vector 'x'. */
int bvar_weight(bvar_t x);

/* Returns first n1 bits of x */
bvar_t bvar_first_bits(bvar_t x, int n1);

/* Returns last n2 bits of x as bvar value */
bvar_t bvar_last_bits(bvar_t x, int n, int n1);

#endif
