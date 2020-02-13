# include "vmcsolver.h"

int main()
{
    VMCSolver *test = new VMCSolver(3, 10, 100);
    test->simulateSystem();
    return 0;
}
