
#include <stdio.h>
#include <stdlib.h>
#include "rand.h"
#include "bjorklund.h"
#include "test_bjorklund.h"

/* Test if the Bjorklund et al.'s parity algorithm works on 'T' random systems of 'm' equations in 'n' unknowns. */
void test_bjorklund_parity(int m, int n, int T)
{
    int c, t;
    for (c = 0, t = 0; t < T; t++)
    {
        qsyst_t *qsyst;

        qsyst = qsyst_new_random(m, n);

        printf("System #%d\n", t + 1);

        double lambda = 0.1967;
        int number_of_iterations = 0;

        int num_sol = qsyst_count(qsyst);

        printf("Number of solutions = %d (parity = %d)\n", num_sol, num_sol % 2);

        int parity = bjorklund_parity(qsyst, lambda, number_of_iterations);

        printf("Randomized algorithm says: parity = %d\n", parity);

        if (parity == num_sol % 2)
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

void test_bjorklund()
{
    int m, n, T = 5;

    printf("Test Bjorklund et al.'s algorithm:\n\n");

    /* Test m = n systems */
    m = 14; n = 14;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_bjorklund_parity(m, n, T);
    printf("\n");

    /* Test m > n systems */
    m = 14; n = 13;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_bjorklund_parity(m, n, T);
    printf("\n");

    /* Test m < n systems */
    m = 13; n = 14;
    printf("Test random %d equations, %d variables systems\n\n", m, n);
    test_bjorklund_parity(m, n, T);
    printf("\n\n");
}
