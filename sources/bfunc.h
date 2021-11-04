#ifndef BFUNC_H
#define BFUNC_H

#include "bvar.h"

/* Boolean function type. It is stored either as:
 * 
 * - A truth table, so for example 
 * 
 * bfunc_get(b0010) is equal to f(0,0,1,0)
 * 
 * - Or a sum of monomials, so for example
 * 
 * b_func_get(b1011) is equal to the coefficient of x_4 x_2 x_1
 * 
 */
typedef struct bfunc_s bfunc_t;

/* Create a new Boolean function of 'n' variables. */
bfunc_t *bfunc_new(int n);

/* Destroy 'bfunc'. */
void bfunc_free(bfunc_t *bfunc);

/* Return a copy of 'bfunc'. */
bfunc_t *bfunc_copy(bfunc_t *bfunc);

/* Get the value or monomial coefficient of the Boolean function at 'x'. */
bool bfunc_get(bfunc_t *bfunc, bvar_t x);

/* Set the value or monomial coefficient of the Boolean function at 'x'. */
void bfunc_set(bfunc_t *bfunc, bvar_t x, bool y);

/* Add 'y' to the value or monomial coefficient of the Boolean function at 'x'. */
void bfunc_add(bfunc_t *bfunc, bvar_t x, bool y);

/* Print the values or monomial coefficients of the Boolean function.
 * (For testing/debugging.) */
void bfunc_print(bfunc_t *bfunc);

/* Replace 'bfunc1' by 'b_func1 & bfunc2'.*/
void bfunc_and(bfunc_t *bfunc1, bfunc_t *bfunc2);

/* Compute the zeta transform of 'bfunc' using Yate's algorithm and
 * overwriting 'bfunc' with the result. */
void bfunc_zeta_transform(bfunc_t *bfunc);

/* Compute the zeta transform of 'bfunc' restricted to words of
 * Hamming weight at most 'd' using Yate's algorithm and 
 * overwriting 'bfunc' with the result. */
void bfunc_restricted_zeta_transform(bfunc_t *bfunc, int d);

/* Compute the zeta transform of 'bfunc' restricted to words of
 * the form (w1, w2), where w1 has length n - n1 and
 * Hamming weight at most 'w' while w2 is a word of length n1, using Yate's algorithm and
 * overwriting 'bfunc' with the result. */
void bfunc_partially_restricted_zeta_transform(bfunc_t *bfunc, int n1, int w);

#endif
