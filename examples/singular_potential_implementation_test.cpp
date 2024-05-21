#include "singular_potential_implementation.hpp"

#include <iostream>

#include "../eigen/Eigen/Dense"

constexpr ball_majumdar_singular_potential::NematicDimension
    dim = ball_majumdar_singular_potential::NematicDimension::full_3D;

using vec = Eigen::Vector<double, ball_majumdar_singular_potential::QTensorShape<dim>::n_degrees_of_freedom>;

int main()
{
    namespace bmsp = ball_majumdar_singular_potential;

    constexpr auto dim = bmsp::NematicDimension::full_3D;
    
    const unsigned int lebedev_order = 590;
    const double damping_parameter = 1.0;
    const double tolerance = 1e-9;
    const unsigned int maximum_iterations = 20;

    bmsp::SingularPotential<dim> singular_potential(lebedev_order,
                                                    damping_parameter,
                                                    tolerance,
                                                    maximum_iterations);

    vec Q = {1.0 / 3.0, -1.0 / 6.0, 0.0, 0.0, 0.0};

    singular_potential.invert_Q(Q);
    auto Lambda = singular_potential.return_Lambda();

    std::cout << Lambda << std::endl;

    return 0;
}
