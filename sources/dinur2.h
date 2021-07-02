/*
 * Interface to the implementation of the algorithms of:
 *
 * Dinur,
 * Cryptanalytic Applications of the Polynomial Method for Solving Multivariate Equation Systems over GF(2)
 * Proceedings Advances in Cryptology - EUROCRYPT 2021.
 *
 */
 
#ifndef DINUR2_H
#define DINUR2_H

#include "qsyst.h"

/* Return True and update 'solution' with up to number_of_solutions solutions of 'qsyst' is found,
 * False otherwise. */
int dinur2_solve(qsyst_t *qsyst, int n1, bvar_t *solution, int limit_number_of_solutions);

#endif
