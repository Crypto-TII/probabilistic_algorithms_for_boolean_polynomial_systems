#ifndef RBFUNC_H
#define RBFUNC_H

#include "bfunc.h"
#include "qsyst.h"

/* Restricted Boolean function type.
 * (Boolean function whose argument has bounded Hamming weight.)
 * 
 *  It is stored either as:
 *
 * - A truth table, so for example
 *
 * rbfunc_get(b0010) is equal to f(0,0,1,0)
 *
 * - Or a sum of monomials, so for example
 *
 * rb_func_get(b1011) is equal to the coefficient of x_4 x_2 x_1
 *
 */

typedef struct rbfunc_s rbfunc_t;

/* Create a new restricted Boolean function of 'n' variables and bounded weight 'w'. */
rbfunc_t *rbfunc_new(int n, int w);

/* Create a new restricted Boolean function representing the polynomial 
 * F = prod_{i=1}^m (p_i + 1), where p_1,...,p_m are the equations of 'qsyst'.
 * This polynomial has the property that F(x) = 1 if and only if 'x' is a solution of 'qsyst'. */
rbfunc_t *rbfunc_new_characteristic(qsyst_t *qsyst);

/* Destroy 'rbfunc'. */
void rbfunc_free(rbfunc_t *rbfunc);

/* Get the value or monomial coefficient of the restricted boolean function at 'x'. */
bool rbfunc_get(rbfunc_t *rbfunc, bvar_t x);

/* Set the value or monomial coefficient of the restricted boolean function at 'x'. */
void rbfunc_set(rbfunc_t *rbfunc, bvar_t x, bool y);

/* Return a copy of 'rbfunc'. */
rbfunc_t *rbfunc_copy(rbfunc_t *rbfunc);

/* Replace 'rbfunc1' with 'rbfunc1 + rbfunc2'. */
void rbfunc_add(rbfunc_t *rbfunc1, rbfunc_t *rbfunc2);

/* Print the values or monomial coefficients of the restricted boolean function.
 * (For testing/debugging.) */
void rbfunc_print(rbfunc_t *rbfunc);

/* replace the value in the position x_index = bvar_get_index(x, rbfunc->n)
 * of rbfunc  x_index xor y . */
void rbfunc_xor(rbfunc_t *rbfunc, bvar_t x, bool y);

/* Compute the zeta transform of in the set W_{w}^{n} 'rbfunc' using Yate's algorithm and
 * overwriting 'rbfunc' with the result. */
void rbfunc_zeta_transform(rbfunc_t *rbfunc);

/* Return a copy of 'rbfunc' as bfunc_t. */
bfunc_t *rbfunc_to_bfunc(rbfunc_t *rbfunc);

#endif
