
#include "consistency.h"
#include "polynomial_solver.h"

bool polynomial_solver(qsyst_t *qsyst, bvar_t *solution, int algorithm_to_use, double e, int number_of_iterations, int up_to_one_solution)
{
    *solution = 0;

    /* Check if 'qsyst' has a zero solution. */
    if (qsyst_are_all_zero(qsyst))
    {
        return 1;
    }
    
    qsyst_t *qsyst1 = qsyst_copy(qsyst);

    int n = qsyst_num_vars(qsyst);
    
    int i;
    for (i = 0; i < n; i++) 
    {
        qsyst_t *qsyst2 = qsyst_copy(qsyst1);
        qsyst_subs(qsyst2, (bvar_t)1, 1);
    
        if (consistency(qsyst2, algorithm_to_use, e, number_of_iterations, up_to_one_solution))
        {
            *solution |= ((bvar_t)1 << (n - i - 1));
            qsyst_subs(qsyst1, (bvar_t)1, 1);
        }
        else
        {
            qsyst_subs(qsyst1, (bvar_t)0, 1);
        }
        
        qsyst_free(qsyst2);
    }
    
    qsyst_free(qsyst1);
    
    if (*solution != 0)
    {
        return 1;
    }

    return 0;
}

