# include <armadillo>
# include <iostream>

using namespace arma;

int main()
{
    mat A(5,5, fill::randu);

    mat L,U,P;

    lu(L,U,P,A);

    L.print("L = ");

}