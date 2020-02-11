#include <iostream>
#include <armadillo>
#include <typeinfo>
#include <math.h> 

using namespace arma;
using namespace std; 

double WaveFunction(double alpha, double beta, arma::mat r, int n_particles, int num_dims, double a);
double CorrelationWaveFunction(arma::subview_row<double> r_i, arma::subview_row<double> r_j, double a, int num_dims);
double LocalEnergy(arma::mat r, double alpha, double beta);


int main(int nargs, char * args[])
{
    if(nargs < 3 || nargs > 3)
    {
        cout << "The file takes two input argument (Number of dimensions, Number of particles)" << endl;
        exit(2);
    }
    int num_particles = atoi(args[2]);
    int num_dimensions = atoi(args[1]);

    mat A(num_particles, num_dimensions, fill::randu);
    // A.print("A: ");

    double alpha = 2.;
    double beta = 1.;
    double a = 0.;

    cout << "WaveFunction value... = " << WaveFunction(alpha, beta, A, num_particles, num_dimensions, a) << endl;
}

/**
 * Function for calculating the value of the wavefunction for a number of given points
 * for given alpha and beta values. 
 * 
 * @param 
 * @returns
 */
double WaveFunction(double alpha, double beta, arma::mat r, int n_particles, int num_dims, double a)
{
    double exp_argument = 0;

    //Calculating the argument of the exponential part of the 
    //wavefunction. 
    for(int i = 0; i < n_particles; i++)
    {

        double temp = 0;

        // Calculating the argument for each particle
        for(int dim = 0; dim < num_dims; dim++)
        {
            if(dim == 2)
            {
                temp += beta*r(i,dim)*r(i,dim);
            }
            else 
            {
                temp += r(i,dim)*r(i,dim);
            }
        }

        exp_argument += temp;
    }

    //Takes the exponential of the argument. 
    double g_product = exp(-alpha*exp_argument);

    //Calculates f (correlation wavefunction) in the wavefunction. 
    double test = CorrelationWaveFunction(r.row(2), r.row(3), num_dims, a);

    double corr_func_product = CorrelationWaveFunction(r.row(0), r.row(1), num_dims, a);

    cout << corr_func_product << endl;

    for (int k = 2; k < n_particles; k++)
    {
        for (int j = 0; j < k; j++)
        {
            corr_func_product *= CorrelationWaveFunction(r.row(j), r.row(k), num_dims, a);
            cout << corr_func_product << endl;
        }
    }

    cout << g_product << endl;
    cout << corr_func_product << endl;

    return g_product*corr_func_product;
}

/**
 * Function for calculating the correlation wavefunction
 * @param
 * @returns
 */
double CorrelationWaveFunction(arma::subview_row<double> r_i, arma::subview_row<double> r_j, double a, int num_dims)
{
    double distance_between_points = 0;

    for(int i = 0; i < num_dims; i++)
    {
        distance_between_points += (r_i(i) - r_j(i))*(r_i(i) - r_j(i));
    }

    distance_between_points = sqrt(distance_between_points);

    double return_value;

    if(abs(distance_between_points) < a)
    {
        return_value = 0;
    }
    else
    {
        return_value = 1 - a/distance_between_points;
    }
    return return_value;
}

/**
 *  Function for calculating the local energy of the system.
 * @param 
 * @returns
 */  
double LocalEnergy(arma::mat r, double alpha, double beta)
{

}