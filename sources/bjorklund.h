/* 
 * Interface to the implementation of the algorithms of
 * 
 * Bjorklund, Kaski, and Williams,
 * Solving systems of polynomial equations over GF(2) by a parity-counting self-reduction,
 * 46th International Colloquium on Automata, Languages, and Programming (ICALP 2019)
 * Article No. 26; pp. 26:1-26:13.
 * 
 */
 
#ifndef BJORKLUND_H
#define BJORKLUND_H

#include "qsyst.h"

/* Return the parity of the number of solutions of 'qsyst'. */
bool bjorklund_parity(qsyst_t* qsyst, double lambda, int number_of_iterations);

#endif
