#include "ball_majumdar_singular_potential.hpp"

#include <iostream>

#include "../eigen/Eigen/Dense"
namespace bmsp = ball_majumdar_singular_potential;

constexpr auto dim = bmsp::NematicDimension::full_3D;

constexpr auto dim2 = bmsp::NematicDimension::quasi_2D;

using vec = Eigen::Vector<double, bmsp::QTensorShape<dim>::n_degrees_of_freedom>;
using vec2 = Eigen::Vector<double, bmsp::QTensorShape<dim2>::n_degrees_of_freedom>;

int main()
{
    const unsigned int lebedev_order = 590;
    const double damping_parameter = 1.0;
    const double tolerance = 1e-9;
    const unsigned int maximum_iterations = 100;

    bmsp::SingularPotential<dim> singular_potential(lebedev_order,
                                                    damping_parameter,
                                                    tolerance,
                                                    maximum_iterations);

    vec Q = {1.0 / 3.0, -1.0 / 6.0, 0.0, 0.0, 0.0};

    singular_potential.invert_Q(Q);
    auto Lambda = singular_potential.return_Lambda();

    bmsp::SingularPotential<dim2> singular_potential2(lebedev_order,
                                                      damping_parameter,
                                                      tolerance,
                                                      maximum_iterations);

    vec2 Q2 = {1.0 / 3.0, -1.0 / 6.0, 0.0};

    singular_potential2.invert_Q(Q2);
    auto Lambda2 = singular_potential.return_Lambda();

    std::cout << Lambda << std::endl;
    std::cout << "\n" << std::endl;
    std::cout << Lambda2 << std::endl;

    return 0;
}

