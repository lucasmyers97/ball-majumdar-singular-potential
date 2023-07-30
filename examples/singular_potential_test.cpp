#include "singular_potential.hpp"
#include "singular_potential.inl"
#include "Q_tensor_shape.hpp"

#include <iostream>

int main()
{
    namespace bmsp = ball_majumdar_singular_potential;

    constexpr auto dim = bmsp::NematicDimension::quasi_2D;
    
    const unsigned int lebedev_order = 590;
    const double damping_parameter = 1.0;
    const double tolerance = 1e-9;
    const unsigned int maximum_iterations = 20;

    bmsp::SingularPotential<dim> singular_potential(lebedev_order,
                                                    damping_parameter,
                                                    tolerance,
                                                    maximum_iterations);

    return 0;
}
