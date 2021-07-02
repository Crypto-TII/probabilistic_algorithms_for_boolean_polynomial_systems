#ifndef POLYNOMIAL_SOLVER_H
#define POLYNOMIAL_SOLVER_H

#include "bvar.h"
#include "qsyst.h"

/* If 'qsyst' has a solution return True and set 'solution',
 * otherwise return False.
 * 
 * algorithm_to_use = 0 uses lokshtanov.c to check the consistency
 * algorithm_to_use = 1 uses bjorklund.c to check the consistency
 * algorithm_to_use = 2 uses dinur.c to check the consistency 
 * 
 * 1 - e is the probability that the isolation algorithm works. */
bool polynomial_solver(qsyst_t *qsyst, bvar_t *solution, int algorithm_to_use, double e, int number_of_iterations, int up_to_one_solution);

#endif
