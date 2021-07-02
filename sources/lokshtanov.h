/* 
 * Interface to the implementation of the algorithms of:
 * 
 * Lokshtanov, Paturi, Tamaki, Williams, and Yu,
 * Beating brute force for systems of polynomial equations over finite fields,
 * Proceedings of the 2017 Annual ACM-SIAM Symposium on Discrete Algorithms
 * 
 */

#ifndef LOKSHTANOV_H
#define LOKSHTANOV_H

#include "qsyst.h"

/* Return True if 'qsyst' has a solution, and False otherwise. */
bool lokshtanov_consistency(qsyst_t *qsyst, double delta, int number_of_iterations);

#endif
