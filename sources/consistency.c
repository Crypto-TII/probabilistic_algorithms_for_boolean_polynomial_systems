
#include <math.h>
#include "consistency.h"
#include "lokshtanov.h"
#include "bjorklund.h"
#include "dinur.h"

int consistency(qsyst_t *qsyst, int algorithm_to_use, double e, int number_of_iterations, int up_to_one_solution)
{
    int n = qsyst_num_vars(qsyst);

    if (n == 0)
    {
        return qsyst_are_all_zero(qsyst);
    }
    
    if (algorithm_to_use == 0)
    {
        double delta = 0.1235;
    
        return lokshtanov_consistency(qsyst, delta, number_of_iterations);
    }
    
    int a = floor(8 * log(1 / e));
    
    if (algorithm_to_use == 1)
    {
        double lambda = 0.1967;

        if (up_to_one_solution)
        {
            if (bjorklund_parity(qsyst, lambda, number_of_iterations)){
                return 1;
            }
            else{
                return 0;
            }
        }
        else{
            int k;
            for (k = 0; k < n; k++)
            {
                int i;
                for (i = 0; i < a; i++)
                {
                    qsyst_t *qsyst_isolated = qsyst_join_linear_polynomials(qsyst, k + 2);

                    if (bjorklund_parity(qsyst_isolated, lambda, number_of_iterations))
                    {
                        qsyst_free(qsyst_isolated);
                        return 1;
                    }

                    qsyst_free(qsyst_isolated);
                }
            }
        }
        return 0;
    }
    
    if (algorithm_to_use == 2)
    {
        float kappa, lambda;
        kappa = 0.31;
        lambda = 0.2;
        if (up_to_one_solution)
        {
            if (dinur_parity(qsyst,  kappa, lambda, number_of_iterations)){
                return 1;
            }
            else{
                return 0;
            }
        }
        else{
            int k;
            for (k = 0; k < n; k++) {
                int i;
                for (i = 1; i < a; i++) {
                    qsyst_t *qsyst_isolated = qsyst_join_linear_polynomials(qsyst, k + 2);

                    if (dinur_parity(qsyst_isolated, kappa, lambda, number_of_iterations)) {
                        qsyst_free(qsyst_isolated);
                        return 1;
                    }

                    qsyst_free(qsyst_isolated);
                }
            }
        }
        return 0;
    }
    
    /* This is just to avoid: 
     * "warning: control reaches end of non-void function". */
    return 0;
}
