#include "quadrature_points.hpp"
#include "singular_potential.hpp"

#include "Q_tensor_shape.hpp"

#include <lebedev_quadrature.hpp>

#include <stdexcept>

namespace ball_majumdar_singular_potential
{

template <NematicDimension dim>
SingularPotential<dim>::
SingularPotential(const unsigned int lebedev_order,
                  const double damping_parameter,
                  const double tolerance,
                  const unsigned int maximum_iterations)
    : damping_parameter(damping_parameter)
    , tolerance(tolerance)
    , maximum_iterations(maximum_iterations)
    , quadrature_points( lebedev::get_order_enum(lebedev_order) )
{
    if (damping_parameter > 1)
        throw std::invalid_argument("Singular potential damping parameter must be <= 1");
}

} // ball_majumdar_singular_potential
