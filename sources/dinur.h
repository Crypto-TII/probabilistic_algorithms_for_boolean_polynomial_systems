/* 
 * Interface to the implementation of the algorithms of:
 * 
 * Dinur,
 * Improved algorithms for solving polynomial systems over GF(2) by multiple parity-counting,
 * Proceedings of the 2021 ACM-SIAM Symposium on Discrete Algorithms (SODA)
 * 
 */

#ifndef DINUR_H
#define DINUR_H

#include "bfunc.h"
#include "qsyst.h"

/* Return the parity of the number of solutions of 'qsyst'. */
bool dinur_parity(qsyst_t *qsyst, double kappa, double lambda, int number_of_iterations);

/* Return the "multiparity" of 'qsyst'. */
bfunc_t *dinur_multiparity(qsyst_t *qsyst, int n1, int w, double lambda, int number_of_iterations);

/* Return the "multiparity" of 'qsyst' computed nonrecursively. */
bfunc_t *dinur_bruteforce_multiparity(qsyst_t * qsyst, int n1, int w);

#endif
