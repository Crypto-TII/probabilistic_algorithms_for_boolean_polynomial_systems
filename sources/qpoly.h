#ifndef QPOLY_H
#define QPOLY_H

#include "bvar.h"
#include "bpoly.h"

/* Type for a quadratic (Boolean) polynomial. */
typedef struct qpoly_s qpoly_t;

/* Create a zero quadratic polynomial in 'n' variables. */
qpoly_t *qpoly_new_zero(int n);

/* Create a random linear polynomial in 'n' variables. */
qpoly_t *qpoly_new_random_linear(int n);

/* Create a random quadratic polynomial in 'n' variables. */
qpoly_t *qpoly_new_random(int n);

/* Return a copy of 'qpoly'. */
qpoly_t *qpoly_copy(qpoly_t *qpoly);

/* Destroy 'qpoly'. */
void qpoly_free(qpoly_t *qpoly);

/* Print expression of 'qpoly'.
 * (For testing/debugging.) */
void qpoly_print(qpoly_t *qpoly);

/* Evaluate 'qpoly' at 'x'. */
bool qpoly_eval(qpoly_t *qpoly, bvar_t);

/* Add the polynomial 'qpoly2' to 'qpoly1'. */
void qpoly_add(qpoly_t *qpoly1, qpoly_t *qpoly2);

/* In 'qpoly' substitute the last 'k' variables with 'z'. */
void qpoly_subs(qpoly_t *qpoly, bvar_t z, int k);

/* Return a copy of 'qpoly' as a 'bpoly_t' polynomial of 'n' variables. */
bpoly_t *qpoly_to_bpoly(qpoly_t *qpoly);

/* Return a copy of 'qpoly' as a 'bfunc_t' polynomial of 'n' variables. */
bfunc_t *qpoly_to_bfunc(qpoly_t *qpoly);

/* Return the constant term of 'qpoly'. */
bool qpoly_costant_term(qpoly_t *qpoly);

#endif
