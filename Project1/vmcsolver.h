#include <armadillo>

using namespace arma; 

/**
 * Class for doing calculations using Variational Monte Carlo on a system of particles. 
 */
class VMCSolver
{
    public:
        VMCSolver(int num_dims, int num_particles, int num_mc_iterations);
        void simulateSystem();

    
    private:
        double calculateLocalEnergyNoInteractions();
        double calculateWaveFunction(arma::mat r);
        double correlationWaveFunction(arma::subview_row<double> r_i, arma::subview_row<double> r_j);

        double internalPotential(double a, arma::subview_row<double> r_i, arma::subview_row<double> r_j);
        double externalPotentialElliptical(mat r);

        double distanceBetweenPoints(arma::subview_row<double> r_i, arma::subview_row<double> r_j);

        void initializeSystemRandom(int num_particles, int num_dimensions);
        void updateSystemConfiguration(arma::mat new_system_configuration);
        arma::mat getNewSystemConfiguration(double step_length);
        arma::mat getSystemConfiguration();

        int _num_dimensions;
        int _num_particles;

        double _alpha;
        double _beta;
        double _a;

        int _num_mc_iterations;

        double _omega_ho = 0.1;
        double _omega_z = 0.1;

        double _mass = 0.000001;

        double _h_bar = 0.00000003;

        arma::mat _system_configuration;  //Matrix of particle positions. Dimensions: (num_particles, num_dimensions)

};

VMCSolver::VMCSolver(int num_dims, int num_particles, int num_mc_iterations) :
    _alpha(0.5), _beta(0.01), _a(0)
{
    _num_dimensions = num_dims;
    _num_particles = num_particles;

    _num_mc_iterations = num_mc_iterations;
}

void VMCSolver::simulateSystem()
{
    initializeSystemRandom(_num_particles, _num_dimensions);
    _system_configuration.print("Before");
    mat A = getNewSystemConfiguration(0.5);
    updateSystemConfiguration(A);
    _system_configuration.print("After");
    cout <<"Wave function value: " << calculateWaveFunction(_system_configuration) << endl;

}

//Initializes the system with given number of particles and dimensions. 
void VMCSolver::initializeSystemRandom(int num_particles, int num_dimensions)
{
    _system_configuration = mat(_num_particles, _num_dimensions, fill::randu);
}

//Returns the current system configuration. 
arma::mat VMCSolver::getSystemConfiguration()
{
    return _system_configuration;
}

//Updates the system system configuration if accepted in the Metropolis algorithm. 
void VMCSolver::updateSystemConfiguration(arma::mat new_system_configuration)
{
    _system_configuration = new_system_configuration; 
}

//Function for returing a possible new configuration. Does NOT set the current configuration to the new one. 
arma::mat VMCSolver::getNewSystemConfiguration(double step_length)
{
    return _system_configuration + step_length*(mat(_num_particles, _num_dimensions, fill::randu) - 0.5);
}

double VMCSolver::calculateWaveFunction(arma::mat r)
{
    double exp_argument = 0;

    //Calculating the argument of the exponential part of the 
    //wavefunction. 
    for(int i = 0; i < _num_particles; i++)
    {

        double temp = 0;

        // Calculating the argument for each particle
        for(int dim = 0; dim < _num_dimensions; dim++)
        {
            if(dim == 2)
            {
                temp += _beta*r(i,dim)*r(i,dim);
            }
            else 
            {
                temp += r(i,dim)*r(i,dim);
            }
        }

        exp_argument += temp;
    }

    //Takes the exponential of the argument. 
    double g_product = exp(-_alpha*exp_argument);

    //Calculates f (correlation wavefunction) in the wavefunction. 
    // double test = correlationWaveFunction(r.row(2), r.row(3));

    double corr_func_product = correlationWaveFunction(r.row(0), r.row(1));

    // cout << corr_func_product << endl;

    for (int k = 2; k < _num_particles; k++)
    {
        for (int j = 0; j < k; j++)
        {
            corr_func_product *= correlationWaveFunction(r.row(j), r.row(k));
            // cout << corr_func_product << endl;
        }
    }

    // cout << g_product << endl;
    // cout << corr_func_product << endl;

    return g_product*corr_func_product;
}

double VMCSolver::correlationWaveFunction(arma::subview_row<double> r_i, arma::subview_row<double> r_j)
{
    double distance_between_points = 0;

    for(int i = 0; i < _num_dimensions; i++)
    {
        distance_between_points += (r_i(i) - r_j(i))*(r_i(i) - r_j(i));
    }

    distance_between_points = sqrt(distance_between_points);

    double return_value;

    if(abs(distance_between_points) < _a)
    {
        return_value = 0;
    }
    else
    {
        return_value = 1 - _a/distance_between_points;
    }
    return return_value;
}

double VMCSolver::calculateLocalEnergyNoInteractions()
{
    double prefix = _h_bar*_h_bar/_mass;

    double temp;

    if(_num_dimensions == 2)
    {
        temp = 2 + _beta;
    }
    else
    {
        temp = _num_dimensions + 1;
    }

    double term_1 = prefix*_alpha*temp*_num_particles; //Stores the value of the first term in the laplacian of phi_k. (-2*alpha)*(2 + beta)*phi_k
    double term_2 = 0; // Stores the value of the second term in the laplacian of phi_k. (4*alpha**2)*(x_k^2 + y_k^2 + beta*z_k^2)*phi_k

    for(int particle; particle < _num_particles; particle++)
    {
        double temp_2 = 0;
        for(int dim; dim < _num_dimensions; dim++)
        {
            if(dim == 2)
            {
                temp_2 += _beta*_beta*_system_configuration(particle, dim)*_system_configuration(particle, dim);
            } 
            else 
            {
                temp_2 += _system_configuration(particle, dim)*_system_configuration(particle, dim);
            }
        }
        term_2 += temp_2;
    }
    term_2 = term_2*prefix*(-2*_alpha*_alpha);

    return term_1 + term_2;
}

double VMCSolver::internalPotential(double a, arma::subview_row<double> r_i, arma::subview_row<double> r_j)
{
    double potential_value;
    double distance_between_points = distanceBetweenPoints(r_i, r_j);

    if(distance_between_points <= a )
    {
        potential_value = 1000000000000000;
    }
    else
    {
        potential_value = 0;
    }
    return potential_value; 
}

//Calculates the external potential (elliptical) for all the particles
double VMCSolver::externalPotentialElliptical(mat r)
{
    double ext_pot = 0;

    for(int particle; particle < _num_particles; particle++)
    {
        double temp = 0;

        for(int dimension; dimension < _num_dimensions; dimension++)
        {
            if(dimension == 2)
            {
                temp += _omega_z*_omega_z*r(particle, dimension)*r(particle, dimension);
            }
            else
            {
                temp += _omega_ho*_omega_ho*r(particle, dimension)*r(particle, dimension);
            }
            
            
        }

        ext_pot += 0.5*_mass;
    }
}

double VMCSolver::distanceBetweenPoints(arma::subview_row<double> r_i, arma::subview_row<double> r_j)
{
    double distance_between_points = 0;

    for(int i = 0; i < _num_dimensions; i++)
    {
        distance_between_points += (r_i(i) - r_j(i))*(r_i(i) - r_j(i));
    }

    distance_between_points = sqrt(distance_between_points);
}