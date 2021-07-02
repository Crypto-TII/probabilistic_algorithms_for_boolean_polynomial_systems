
#include <stdlib.h>
#include "bfunc.h"
#include "rbfunc.h"
#include "dinur2.h"


void compute_uvalues(qsyst_t *qsyst, int n1, int w, rbfunc_t **ZV)
{
    int n = qsyst_num_vars(qsyst);
    bvar_t x, x2, y, z, x1 = 0;
    int h = 0;
    do
    {
        for (x2 = 0; x2 < (1 << n1); x2++)
        {
            x = (x1 << n1) | x2;
            if (qsyst_is_solution(qsyst, x))
            {
                z = bvar_first_bits(x, n1);
                y = bvar_last_bits(x, n, n1);
                if (bvar_weight(y) <= w)
                {
                    rbfunc_xor(ZV[0], y, 1);
                }
                int j;
                for (j = 1; j < n1 + 1; j++)
                {
                    if ((z & (1 << (n1 - j))) == 0)
                    {
                        rbfunc_xor(ZV[j], y, 1);
                    }
                }
            }
        }
    }
    while (next_subset(&x1, &h, n - n1, w + 1));
}

void output_potential_solutions(qsyst_t *qsyst, int n1, int w, bfunc_t **current_potential_outputs)
{
    int n = qsyst_num_vars(qsyst);

    rbfunc_t **ZV = malloc((n1 + 1) * sizeof(rbfunc_t*));

    ZV[0] = rbfunc_new(n - n1, w);
    int k;
    for (k = 1; k < n1 + 1; k ++)
    {
        ZV[k] = rbfunc_new(n - n1, w + 1);
    }
    compute_uvalues(qsyst, n1, w, ZV);

    for (k = 0; k < n1 + 1; k ++){
        rbfunc_zeta_transform(ZV[k]);
    }

    bfunc_t **evals;
    evals = malloc((n1 + 1) * sizeof(bfunc_t*));

    for (k = 0; k < n1 + 1; k++)
    {
        evals[k] = rbfunc_to_bfunc(ZV[k]);
        rbfunc_free(ZV[k]);
        bfunc_zeta_transform(evals[k]);
    }
    free(ZV);
    bvar_t y;
    for (y = 0; y < (1 << (n - n1)); y++)
    {
        if (bfunc_get(evals[0], y))
        {
            bfunc_set(current_potential_outputs[0], y, 1);
            int i;
            for (i = 1; i < n1 + 1; i++)
            {
                bfunc_set(current_potential_outputs[i], y, !bfunc_get(evals[i], y));
            }
        }
    }
    int i;
    for (i = 0; i < n1 + 1; i++)
    {
        bfunc_free(evals[i]);
    }
    free(evals);
}

int dinur2_solve(qsyst_t *qsyst, int n1, bvar_t *solution, int limit_number_of_solutions)
{
    int n = qsyst_num_vars(qsyst);
    int w = n1 + 2;
    int N = 5;

    bfunc_t **potential_solutions_list[N];
    if ((n - n1) <= w + 1)
    {
        int c = qsyst_bruteforce(qsyst, solution, limit_number_of_solutions);
        return c;
    }

    int i, k;
    int c = 0;
    for (k = 0; k < N; k++)
    {
        bfunc_t **current_potential_outputs;
        current_potential_outputs = malloc((n1 + 1) * sizeof(bfunc_t*));
        for (i = 0; i < n1 + 1; i++)
        {
            current_potential_outputs[i] = bfunc_new(n - n1);
        }

        qsyst_t *qsyst_lin_com = qsyst_rand_lin_comb(qsyst, n1 + 1);
        output_potential_solutions(qsyst_lin_com, n1, w, current_potential_outputs);
        potential_solutions_list[k] = current_potential_outputs;

        bvar_t y;
        bvar_t y_max = 1 << (n - n1);
        for (y = 0; y < y_max; y++)
        {
            if (bfunc_get(current_potential_outputs[0], y))
            {
                int k1;
                for (k1 = 0; k1 < k; k1++)
                {
                    int k3;
                    int solution_found_before = 1;
                    if (bfunc_get(potential_solutions_list[k1][0], y) != 1)
                    {
                        solution_found_before = 0;
                    }
                    for (k3 = 1; k3 < n1 + 1; k3++)
                    {
                        if (bfunc_get(current_potential_outputs[k3], y) !=
                            bfunc_get(potential_solutions_list[k1][k3], y))
                        {
                            solution_found_before = 0;
                        }
                    }

                    if (solution_found_before)
                    {
                        int k2;
                        bvar_t sol = y << n1;
                        for (k2 = 1; k2 < n1 + 1; k2++)
                        {
                            sol |= bfunc_get(current_potential_outputs[k2], y) << (n1 - k2);
                        }

                        if (qsyst_is_solution(qsyst, sol))
                        {
                            bool sol_was_found_before = 0;
                            int l;
                            for (l = 0; l < c; l ++)
                            {
                                if (sol == solution[l])
                                {
                                    sol_was_found_before = 1;
                                }

                            }
                            if (!sol_was_found_before)
                            {
                                solution[c] = sol;
                                c += 1;
                            }

                            if (c == limit_number_of_solutions)
                            {
                                int i1, i2;
                                for (i1 = 0; i1 < k + 1; i1++)
                                {
                                    for (i2 = 0; i2 < n1 + 1; i2++)
                                    {
                                        bfunc_free(potential_solutions_list[i1][i2]);
                                    }
                                }

                                for (i1 = 0; i1 < k + 1; i1++){
                                    free(potential_solutions_list[i1]);
                                }

                                qsyst_free(qsyst_lin_com);
                                return c;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }
        qsyst_free(qsyst_lin_com);
    }
    int i1, i2;
    for (i1 = 0; i1 < N; i1++)
    {
        for (i2 = 0; i2 < n1 + 1; i2++)
        {
            bfunc_free(potential_solutions_list[i1][i2]);
        }
    }

    for (i1 = 0; i1 < N; i1++)
    {
        free(potential_solutions_list[i1]);
    }

    return c;
}



