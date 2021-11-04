
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rand.h"
#include "bfunc.h"
#include "qpoly.h"
#include "qsyst.h"

struct qsyst_s
{
    int           m; /* Number of equations */
    int           n; /* Number of variables */
    qpoly_t **qpoly; /* Equations */
};

qsyst_t *qsyst_new_random(int m, int n)
{
    qsyst_t *qsyst;
    
    qsyst = malloc(sizeof(qsyst_t));
    qsyst->m = m;
    qsyst->n = n;
    qsyst->qpoly = malloc(m * sizeof(qpoly_t*));
    
    int i;
    for (i = 0; i < m; i++)
    {
        qsyst->qpoly[i] = qpoly_new_random(n);
    }

    return qsyst;    
}

void qsyst_free(qsyst_t *qsyst)
{
    int i;
    for (i = 0; i < qsyst->m; i++)
    {
        qpoly_free(qsyst->qpoly[i]);
    }
    
    free(qsyst->qpoly);
    
    free(qsyst);
}

int qsyst_num_vars(qsyst_t *qsyst)
{
    return qsyst->n;
}

int qsyst_num_eqs(qsyst_t *qsyst)
{
    return qsyst->m;
}

qpoly_t *qsyst_equ(qsyst_t *qsyst, int j)
{
    return qsyst->qpoly[j];
}

void qsyst_print(qsyst_t *qsyst)
{
    int i;
    for (i = 0; i < qsyst->m; i++)
    {
        printf("(%d) : ", i + 1);
        qpoly_print(qsyst->qpoly[i]);
    }
}

bool qsyst_is_solution(qsyst_t *qsyst, bvar_t x)
{
    int i = 0;
    while (i < qsyst->m && qpoly_eval(qsyst->qpoly[i], x) == 0)
    {
        i++;
    }
    
    return i >= qsyst->m;
}

int qsyst_count(qsyst_t *qsyst)
{
    int c;
    bvar_t x;
    for (x = 0, c = 0; x < ((bvar_t)1 << qsyst->n); x++)
    {
        if (qsyst_is_solution(qsyst, x))
        {
            c++;
        }
    }
    
    return c;
}

qsyst_t *qsyst_copy(qsyst_t *qsyst)
{
    qsyst_t *copy;
    copy = malloc(sizeof(qsyst_t));
    copy->m = qsyst->m;
    copy->n = qsyst->n;
    copy->qpoly = malloc(qsyst->m * sizeof(qpoly_t*));
    
    int i;
    for (i = 0; i < qsyst->m; i++)
    {
        copy->qpoly[i] = qpoly_copy(qsyst->qpoly[i]);
    }
    
    return copy;
}

void qsyst_subs(qsyst_t *qsyst, bvar_t z, int k)
{
    int i;
    for (i = 0; i < qsyst->m; i++)
    {
        qpoly_subs(qsyst->qpoly[i], z, k);
    }
    qsyst->n -= k;
}

qsyst_t *qsyst_rand_lin_comb(qsyst_t* qsyst, int r)
{
    qsyst_t* rlc;
    rlc = malloc(sizeof(qsyst_t));
    rlc->m = r;
    rlc->n = qsyst->n;
    rlc->qpoly = malloc(r * sizeof(qpoly_t*));
    
    int i, j;
    for (i = 0; i < r; i++)
    {
        rlc->qpoly[i] = qpoly_new_zero(qsyst->n);
        
        for (j = 0; j < qsyst->m; j++)
        {
            if (rand_bool())
            {
                qpoly_add(rlc->qpoly[i], qsyst->qpoly[j]);
            }
        }
    }

    return rlc;
}

qsyst_t *qsyst_join_linear_polynomials(qsyst_t *qsyst, int r)
{
    qsyst_t *qsyst_output;

    int n, m;
    n = qsyst_num_vars(qsyst);
    m = qsyst_num_eqs(qsyst);

    qsyst_output = malloc(sizeof(qsyst_t));
    qsyst_output->m = m + r;
    qsyst_output->n = n;
    qsyst_output->qpoly = malloc((m + r) * sizeof(qpoly_t*));

    int i;
    for (i = 0; i < m; i++)
    {
        qsyst_output->qpoly[i] = qpoly_copy(qsyst->qpoly[i]);
    }

    for (i = 0; i < r; i++)
    {
        qsyst_output->qpoly[m + i] = qpoly_new_random_linear(n);
    }
    
    return qsyst_output;
}

bool qsyst_are_all_zero(qsyst_t *qsyst)
{
    int i;
    for(i = 0; i < qsyst_num_eqs(qsyst); i++)
    {
        qpoly_t *poly = qsyst->qpoly[i];

        if (qpoly_costant_term(poly) == 1)
        {
            return 0;
        }
    }

    return 1;
}

int qsyst_bruteforce(qsyst_t *qsyst, bvar_t *solution, int limit_number_of_solutions)
{
    int n = qsyst_num_vars(qsyst);
    bvar_t x;
    int c = 0;
    for (x = 0; x < ((bvar_t)1 << n); x++)
    {
        if (qsyst_is_solution(qsyst, x))
        {
            solution[c] = x;
            c++;
            if (c == limit_number_of_solutions)
            {
                return c;
            }
        }
    }
    return c;
}
