#include "test_rbfunc.h"
#include "rbfunc.h"
#include "qsyst.h"
#include <stdio.h>
#include <stdlib.h>

void test_rbfunc()
{
    int num_test_passed = 0;
    //Test 1
    int output = bvar_number_of_values(11,  5);
    int expected = 1024;
    if (output == expected)
    {
        num_test_passed += 1;
    }

    //Test 2
    int expected2 = 0;
    rbfunc_t  *rbfunc = rbfunc_new(11, 5);
    int output2 = rbfunc_get(rbfunc, 120);
    rbfunc_free(rbfunc);
    if (output2 == expected2)
    {
        num_test_passed += 1;
    }

    //Test 3
    int expected3 = 1;
    rbfunc_t  *rbfunc3 = rbfunc_new(17, 5);
    rbfunc_set(rbfunc3, 140, 1);
    rbfunc_set(rbfunc3, 10, 1);
    rbfunc_set(rbfunc3, 200, 1);
    int output3 = rbfunc_get(rbfunc3, 10);
    if (output3 == expected3)
    {
        num_test_passed += 1;
    }

    //Test 4
    int expected4 = 0;
    rbfunc_set(rbfunc3, 200, 0);
    int output4 = rbfunc_get(rbfunc3, 200);
    if (output4 == expected4)
    {
        num_test_passed += 1;
    }

    //Test 5
    rbfunc_t *copy_rbfunc3 = rbfunc_copy(rbfunc3);
    int expected5 = 1;
    int output5 = rbfunc_get(copy_rbfunc3, 10);
    if (output5 == expected5)
    {
        num_test_passed += 1;
    }

    //Test 6
    int expected6 = 0;
    int output6= rbfunc_get(copy_rbfunc3, 200);
    if (output6 == expected6)
    {
        num_test_passed += 1;
    }
    rbfunc_free(copy_rbfunc3);
    rbfunc_free(rbfunc3);

    //Test 7
    int expected7 = 1;
    int n = 9, w = 10;
    rbfunc_t  *rbfunc7 = rbfunc_new(n, w);
    bvar_t x = (1 << (w - 1)) - 5;
    rbfunc_set(rbfunc7, x, 1);
    int output7 = rbfunc_get(rbfunc7, x);
    if (output7 == expected7)
    {
        num_test_passed += 1;
    }
    rbfunc_free(rbfunc7);
    int total_tests = 7;

    if(total_tests == num_test_passed)
    {
        printf("All the test passed!\n");
    }
    else
    {
        printf("SOME TEST FAILED!\n");
    }
}

void test_rbfunc_new_characteristic()
{
    int m; int n;
    
    for (n = 10; n < 15; n++)
    {
        for (m = n - 2; m < n + 3; m++)
        {
            printf("m = %d, n = %d\n", m, n);
            
            qsyst_t *qsyst;
            rbfunc_t *rbfunc;
            
            qsyst = qsyst_new_random(m, n);
            
            rbfunc = rbfunc_new_characteristic(qsyst);
            rbfunc_zeta_transform(rbfunc);
        
            bvar_t x;
            for (x = 0; x < ((bvar_t)1 << n); x++)
            {
                if (qsyst_is_solution(qsyst, x) != rbfunc_get(rbfunc, x))
                {
                    printf("Error!\n");
                    exit(0);
                }
            }
            
            rbfunc_free(rbfunc);
            qsyst_free(qsyst);
        }
    }
}
