
#include <stdlib.h>
#include <stdio.h>
#include "rand.h"
#include "test_lokshtanov.h"
#include "lokshtanov.h"

/* Test if the Lokshtanov et al.'s consistency algorithm works on 'T' random systems of 'm' equations in 'n' unknowns. */
void test_lokshtanov_consistency(int m, int n, int T)
{
    int c, t;
    for (c = 0, t = 0; t < T; t++)
    {
        double delta = 0.1235;
        int number_of_iterations = 0;
        qsyst_t *qsyst;

        qsyst = qsyst_new_random(m, n);

        printf("System #%d\n", t + 1);

        //qsyst_print(qsyst);

        int num_sol = qsyst_count(qsyst);

        printf("Number of solutions = %d\n", num_sol);
        int consistency = lokshtanov_consistency(qsyst, delta, number_of_iterations);

        printf("Randomized algorithm says: consistency = %d\n", consistency);

        if ((num_sol > 0) == consistency)
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

void test_lokshtanov()
{
    int m, n, T = 5;
    
    printf("Test Lokshtanov et al.'s algorithm:\n\n");

    /* Test m = n systems */
    m = 10; n = 10;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_lokshtanov_consistency(m, n, T);
    printf("\n");

    /* Test m > n systems */
    m = 10; n = 0;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_lokshtanov_consistency(m, n, T);
    printf("\n");

    /* Test m < n systems */
    m = 9; n = 10;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_lokshtanov_consistency(m, n, T);
    printf("\n\n");
}
