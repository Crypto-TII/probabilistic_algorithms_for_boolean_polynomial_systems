
#include <stdlib.h>
#include <sys/random.h>
#include "rand.h"

void rand_init(int seed)
{
	srand(seed);
}

bool rand_bool()
{
	return rand() & 1;
}
