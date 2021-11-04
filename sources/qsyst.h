#ifndef QSYST_H
#define QSYST_H

#include "bvar.h"
#include "qpoly.h"

/* Type for a system of quadratic polynomials. */
typedef struct qsyst_s qsyst_t;

/* Create a random quadratic system of 'm' equations in 'n' variables. */
qsyst_t *qsyst_new_random(int m, int n);

/* Destroy 'qsyst'. */
void qsyst_free(qsyst_t *qsyst);

/* Return the number of unknowns of 'qsyst'. */
int qsyst_num_vars(qsyst_t *qsyst);

/* Return the number of equations of 'qsyst'. */
int qsyst_num_eqs(qsyst_t *qsyst);

/* Return the 'j'th equation of 'qsyst'. */
qpoly_t *qsyst_equ(qsyst_t *qsyst, int j);

/* Print the system. */
void qsyst_print(qsyst_t *qsyst);

/* Check if 'x' is a solution of 'qsyst'. */
bool qsyst_is_solution(qsyst_t *qsyst, bvar_t x);

/* Return the number of solutions of 'qsyst' (computed by brute force). */
int qsyst_count(qsyst_t *qsyst);

/* Return a copy of 'qsyst'. */
qsyst_t *qsyst_copy(qsyst_t *qsyst);

/* Substitute the last 'k' variables of 'qsyst' with 'z' */
void qsyst_subs(qsyst_t *qsyst, bvar_t z, int k);

/* Create a system of 'r' random linear combinations of the
 * equations of 'qsyst' */
qsyst_t *qsyst_rand_lin_comb(qsyst_t* qsyst, int r);

/* Return a new system equal to the join of 'qsyst' and 'r' random linear polynomials. */
qsyst_t *qsyst_join_linear_polynomials(qsyst_t *qsyst, int r);

/* Return True if all the constant polynomials in 'qsyst' are zero. */
bool qsyst_are_all_zero(qsyst_t *qsyst);

/* Compute a solution of 'qsyst' by bruteforce.
 * Return True if a solution exists and put it in 'sol',
 * False otherwise. */
int qsyst_bruteforce(qsyst_t *qsyst, bvar_t *solution, int number_of_solutions);

#endif
