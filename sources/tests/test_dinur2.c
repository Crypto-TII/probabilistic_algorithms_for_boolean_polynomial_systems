
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tests/test_dinur2.h"
#include "dinur2.h"
#include "rand.h"
#include "qsyst.h"

/* Test if the Dinur's second algorithm works on 'T' random systems of 'm' equations in 'n' unknowns. */
void test_dinur2_solver(int m, int n, int T)
{
    int c, t;
    for (c = 0, t = 0; t < T; t++)
    {
        qsyst_t *qsyst;

        qsyst = qsyst_new_random(m, n);

        printf("System #%d\n", t + 1);

        int num_sol = qsyst_count(qsyst);

        printf("Number of solutions = %d \n", num_sol);

        int limit_number_of_solutions = num_sol;
        bvar_t sol[limit_number_of_solutions];

        int n1 = fmax(floor(n/.54), 1);
        
        int num_solutions_found = dinur2_solve(qsyst, n1, sol, limit_number_of_solutions);
        
        if (num_solutions_found)
        {
            printf("Solutions found: \n");
            int i;
            for (i = 0; i < num_solutions_found; i++)
            {
                bvar_print(sol[i], n);
            }
        }
        else
        {
            printf("No solution found.\n");
        }
        
        bool right =  (num_solutions_found) || (!num_solutions_found && num_sol == 0);
        int i;
        for (i = 0; i < num_solutions_found; i++)
        {
            if (!qsyst_is_solution(qsyst, sol[i]))
            {
                right = 0;
            }
        }
        
        if (right)
        {
            printf("---> Right!\n");
            c++;
        }
        else
        {
            printf("---> Wrong!\n");
            //exit(0);
        }
        printf("\n");
        
        qsyst_free(qsyst);
    }
    
    printf("Guessed %d out of %d (%.2lf%%)\n", c, T, 100 * c / (double)T);

}

void test_dinur2()
{
    int m, n, T = 5;

    printf("Test Dinur's second algorithm:\n\n");

    /* Test m = n systems */
    m = 19; n = 19;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_dinur2_solver(m, n, T);
    printf("\n");

    /* Test m > n systems */
    m = 19; n = 18;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_dinur2_solver(m, n, T);
    printf("\n");

    /* Test m < n systems */
    m = 18; n = 19;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_dinur2_solver(m, n, T);
    printf("\n\n");
}
