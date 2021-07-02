
#include <stdlib.h>
#include <stdio.h>
#include "polynomial_solver.h"
#include "rand.h"
#include "test_polynomial_solver.h"

/* Test if the polynomial solver algorithm works on 'T' random systems of 'm' equations in 'n' unknowns. */
void test_solver(int m, int n, int T)
{
    int t;

    for (t = 0; t < T; t++)
    {
        qsyst_t *qsyst;

        qsyst = qsyst_new_random(m, n);

        printf("System #%d\n", t + 1);

        int num_sol = qsyst_count(qsyst);

        printf("Number of solutions = %d \n", num_sol);

        const char* algo[] = {"Lokshtanov", "Bjorklund", "Dinur"};
        int a;
        int number_of_iterations = 0;
        int up_to_one_solution = 0;
        for (a = 0; a < 3; a++)
        {
            bvar_t sol;
            bool right;
            
            if (polynomial_solver(qsyst, &sol, a, 1.0 / n, number_of_iterations, up_to_one_solution))
            {
                right = qsyst_is_solution(qsyst, sol);
            }
            else
            {
                right = (num_sol == 0); 
            }
            
            if (right)
            {
                printf("%s is right!\n", algo[a]);
            }
            else
            {
                printf("%s is WRONG!\n", algo[a]);
                //exit(0);
            }
        }
        
        printf("\n");
        qsyst_free(qsyst);
    }

}

void test_polynomial_solver()
{
    int m, n, T = 5;

    printf("Test Polynomial_Solver algorithm:\n\n");

    /* Test m = n systems */
    m = 10; n = 10;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_solver(m, n, T);
    printf("\n");

    /* Test m > n systems */
    m = 10; n = 9;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_solver(m, n, T);
    printf("\n");

    /* Test m < n systems */
    m = 9; n = 10;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_solver(m, n, T);
    printf("\n\n");
}
