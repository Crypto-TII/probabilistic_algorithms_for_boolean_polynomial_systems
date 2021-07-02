#ifndef BPOLY_H
#define BPOLY_H

#include "bfunc.h"

/* Structure for a sparse Boolean polynomial. */
typedef struct bpoly_s bpoly_t;

/* Return a new Boolean polynomial equal to 0. */
bpoly_t *bpoly_new_zero();

/* Return a new Boolean polynomial equal to 1. */
bpoly_t *bpoly_new_one();

/* Create a new Boolean polynomial from the array 'mon' of 'num' monomials. */
bpoly_t *bpoly_new(bvar_t *mon, int num);

/* Create a new random Boolean polynomial of 'n' variables and at most 't' monomials. 
 * (For testing/debugging.) */
bpoly_t *bpoly_new_random(int n, int t);

/* Destroy 'bpoly'. */
void bpoly_free(bpoly_t *bpoly);

/* Return a copy of 'bpoly'. */
bpoly_t *bpoly_copy(bpoly_t *bpoly);

/* Print expression of 'bpoly'.
 * (For testing/degugging.) */
void bpoly_print(bpoly_t *bpoly);

/* Return the degree of 'bpoly' (-1 if 'bpoly' is zero). */
int bpoly_degree(bpoly_t *bpoly);

/* Return the variables appearing in 'bpoly'. */
bvar_t bpoly_vars(bpoly_t *bpoly);

/* Return the value of 'bpoly' at 'x'. */
bool bpoly_eval(bpoly_t *bpoly, bvar_t x);

/* Add the monomial 'x' to the polynomial 'bpoly'. */
void bpoly_add_mon(bpoly_t *bpoly, bvar_t x);

/* Multiply the polynomial 'bpoly' for the monomial 'x'. */
void bpoly_mul_mon(bpoly_t *bpoly, bvar_t x);

/* Replace 'bpoly1 with 'bpoly1 + bpoly2'. */
void bpoly_add(bpoly_t *bpoly1, bpoly_t *bpoly2);

/* Replace 'bpoly1 with 'bpoly1 * bpoly2'. */
void bpoly_mul(bpoly_t *bpoly1, bpoly_t *bpoly2);

/* Replace the variables 'x' of 'bpoly' with the values in 'a'.
 * For example, 'bpoly_assign(bpoly, b1011, b1001)' sets
 * x_1 = 1, x_2 = 0, x_4 = 1. */
void bpoly_assign(bpoly_t *bpoly, bvar_t x, bvar_t a);

/* Return the monomials of 'bpoly' as a Boolean function in 'n' variables. */
bfunc_t *bpoly_monomials(bpoly_t *bpoly, int n);

/* Return the truth table of 'bpoly' as a Boolean function in 'n' variables.
 * The truth table is evaluated using the restricted zeta transform. */
bfunc_t *bpoly_truth_table(bpoly_t *bpoly, int n);

#endif
