#ifndef RAND_H
#define RAND_H

#include <stdint.h>
#include <stdbool.h>

/* Initialize pseudorandom number generator with 'seed'. */
void rand_init(int seed);

/* Return a random Boolean. */
bool rand_bool();

#endif
