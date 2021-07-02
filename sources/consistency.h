#ifndef CONSISTENCY_H
#define CONSISTENCY_H

#include "qsyst.h"

/* Return True if the system 'qsyst' has a solution and False otherwise
 * 
 * algorithm_to_use = 0 uses lokshtanov.c to check the consistency
 * algorithm_to_use = 1 uses bjorklund.c to check the consistency
 * algorithm_to_use = 2 uses dinur.c to check the consistency 
 * 
 * 1 - e is the probability that the isolation algorithm works. */
int consistency(qsyst_t *qsyst, int algorithm_to_use, double e, int number_of_iterations, int up_to_one_solution);

#endif
