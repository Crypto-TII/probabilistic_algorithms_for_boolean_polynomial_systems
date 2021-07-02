
#include <stdlib.h>
#include <stdio.h>
#include "tests/test_dinur.h"
#include "dinur.h"
#include "rand.h"
#include <time.h>
/* Test if the Dinur's parity algorithm works on 'T' random systems of 'm' equations in 'n' unknowns. */
void test_dinur_parity(int m, int n, int T)
{
    int c, t;
    
    double kappa = 0.31;
    double lambda = 0.2;
    int number_of_iterations = 0;
    
    for (c = 0, t = 0; t < T; t++)
    {
        qsyst_t *qsyst;

        qsyst = qsyst_new_random(m, n);

        printf("System #%d\n", t + 1);

        int num_sol = qsyst_count(qsyst);
        printf("Number of solutions = %d (parity = %d)\n", num_sol, num_sol % 2);

        clock_t t;
        t = clock();
        int parity = dinur_parity(qsyst, kappa, lambda, number_of_iterations);
        t = clock() - t;
        double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

        printf("Dinur1 took %f seconds to execute \n", time_taken);

        printf("Randomized algorithm says: parity = %d\n", parity);

        if (parity == (num_sol % 2))
        {
            printf("---> Right!\n");
            c++;
        }
        else
        {
            printf("---> Wrong!\n");
            exit(0);
        }

        printf("\n");

        qsyst_free(qsyst);
    }

    printf("Guessed %d out of %d (%.2lf%%)\n", c, T, 100 * c / (double)T);
}

void test_dinur()
{
    int m, n, T = 5;
    
    printf("Test Dinur's algorithm:\n\n");

    /* Test m = n systems */
    m = 14; n = 14;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_dinur_parity(m, n, T);
    printf("\n");

    /* Test m > n systems */
    m = 14; n = 13;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_dinur_parity(m, n, T);
    printf("\n");

    /* Test m < n systems */
    m = 13; n = 14;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_dinur_parity(m, n, T);
    printf("\n\n");
}
