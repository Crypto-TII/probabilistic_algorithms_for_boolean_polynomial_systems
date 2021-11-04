
#include "rand.h"
#include "tests/test_lokshtanov.h"
#include "tests/test_bjorklund.h"
#include "tests/test_dinur.h"
#include "tests/test_dinur2.h"
#include "tests/test_polynomial_solver.h"

int main(int argc, char *argv[])
{
    rand_init(1234321);
    
    test_lokshtanov();
    
    test_bjorklund();
    
    test_dinur();
    
    test_dinur2();
    
    test_polynomial_solver();

    return 0;
}
