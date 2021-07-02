
#include <stdio.h>
#include <stdlib.h>
#include "bpoly.h"

/* Test sums and products of random polynomials of 'n' variables and at most 't' monomials. */
void test_operations(int n, int t)
{
	printf("Testing sums and products...\n");
	
	int k;
	for (k = 0; k < 100; k++)
	{
		bpoly_t *bpoly1, *bpoly2, *sum, *prod;
		
		bpoly1 = bpoly_new_random(n, t);
		bpoly2 = bpoly_new_random(n, t);
		
		sum = bpoly_copy(bpoly1);
		bpoly_add(sum, bpoly2);

		prod = bpoly_copy(bpoly1);
		bpoly_mul(prod, bpoly2);
		
		//printf("P1      = "); bpoly_print(bpoly1);
		//printf("P2      = "); bpoly_print(bpoly2);
		//printf("P1 + P2 = "); bpoly_print(sum);
		//printf("P1 * P2 = "); bpoly_print(prod);
		//printf("\n");
		
		bvar_t x;
		for (x = 0; x < ((bvar_t)1 << n); x++)
		{
			bool v1, v2, s, p;
			
			v1 = bpoly_eval(bpoly1, x);
			v2 = bpoly_eval(bpoly2, x);
			s  = bpoly_eval(sum, x);
			p  = bpoly_eval(prod, x);
			
			if (s != (v1 ^ v2) || p != (v1 & v2))
			{
				printf("ERROR!\n");
				return;
			}
		}
		
		bpoly_free(bpoly1);
		bpoly_free(bpoly2);
		bpoly_free(sum);
		bpoly_free(prod);
	}
	
	printf("Everything is OK!\n");
}

/* Test assignment of random polynomials of 'n' variables and at most 't' monomials. */
void test_assignment(int n, int t)
{
	printf("Testing assignment...\n");
	
	int k;
	for (k = 0; k < 100; k++)
	{
		bpoly_t *bpoly = bpoly_new_random(n, t);
		
		//bpoly_print(bpoly);

		bvar_t v, a;
			
		v = rand() & (((bvar_t)1 << n) - 1);
		a = rand() & (((bvar_t)1 << n) - 1);
		
		bpoly_t *spec = bpoly_copy(bpoly);
		bpoly_assign(spec, v, a);
			
		//bpoly_print(spec);
						
		bvar_t x;
		for (x = 0; x < ((bvar_t)1 << n); x++)
		{
			if ((x & v) == a 
				&& bpoly_eval(bpoly, x) != bpoly_eval(spec, x))
			{
				printf("ERROR!\n");
				return;				
			}
		}
		
		bpoly_free(bpoly);
		bpoly_free(spec);
	}

	printf("Everything is OK!\n");
}

/* Test fast evaluation of a random polynomial in 'n' variables and at most 't' monomials. */
void test_fast_evaluation(int n, int t)
{
	printf("Testing fast evaluation...\n");
	
	int k;
	for (k = 0; k < 100; k++)
	{
		bpoly_t *bpoly;
		bfunc_t *bfunc;
		
		bpoly = bpoly_new_random(n, t);
		
		//bpoly_print(bpoly);
		
		bfunc = bpoly_monomials(bpoly, n);
		
		//bfunc_print(bfunc);
		
		//bfunc_zeta_transform(bfunc);
		bfunc_zeta_transform(bfunc);
	
		bvar_t x;
		for (x = 0; x < ((bvar_t)1 << n); x++)
		{			
			if (bpoly_eval(bpoly, x) != bfunc_get(bfunc, x))
			{
				printf("ERROR!\n");
				return;
			}
		}
		
		bpoly_free(bpoly);
	}

	printf("Everything is OK!\n");
}

void test_bpoly()
{	
	test_operations(5, 12);
	
	test_assignment(5, 12);

	test_fast_evaluation(5, 12);
}
