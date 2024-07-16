#include "preprocessor.hpp"
#include "singular_potential.hpp"

#include "Q_tensor_shape.hpp"

#include "../eigen/Eigen/Dense"
#include "../lebedev-quadrature/src/lebedev_quadrature.hpp"

#include <cmath>
#include <stdexcept>

namespace ball_majumdar_singular_potential
{

template <>
SINGULAR_POTENTIAL_EXTERNAL_LINKAGE
SingularPotential<NematicDimension::full_3D>::
SingularPotential(const unsigned int lebedev_order,
                  const double damping_parameter,
                  const double tolerance,
                  const unsigned int maximum_iterations)
    : damping_parameter(damping_parameter)
    , tolerance(tolerance)
    , maximum_iterations(maximum_iterations)
    , quadrature_points( lebedev::get_order_enum(lebedev_order) )
    , x(quadrature_points.get_x())
    , y(quadrature_points.get_y())
    , z(quadrature_points.get_z())
    , w(quadrature_points.get_weights())
{
    if (damping_parameter > 1)
        throw std::invalid_argument("Singular potential damping parameter must be <= 1");

    delta_vec = {1, 1, 0, 0, 0};
}



template <>
SINGULAR_POTENTIAL_EXTERNAL_LINKAGE
SingularPotential<NematicDimension::quasi_2D>::
SingularPotential(const unsigned int lebedev_order,
                  const double damping_parameter,
                  const double tolerance,
                  const unsigned int maximum_iterations)
    : damping_parameter(damping_parameter)
    , tolerance(tolerance)
    , maximum_iterations(maximum_iterations)
    , quadrature_points( lebedev::get_order_enum(lebedev_order) )
    , x(quadrature_points.get_x())
    , y(quadrature_points.get_y())
    , z(quadrature_points.get_z())
    , w(quadrature_points.get_weights())
{
    if (damping_parameter > 1)
        throw std::invalid_argument("Singular potential damping parameter must be <= 1");

    delta_vec = {1, 1, 0};
}



template <NematicDimension dim>
SINGULAR_POTENTIAL_EXTERNAL_LINKAGE
Eigen::Vector<double, QTensorShape<dim>::n_degrees_of_freedom>
SingularPotential<dim>::
return_Lambda() const
{
    return Lambda;
}



template <NematicDimension dim>
SINGULAR_POTENTIAL_EXTERNAL_LINKAGE
Eigen::Matrix<double, 
              QTensorShape<dim>::n_degrees_of_freedom,
              QTensorShape<dim>::n_degrees_of_freedom>
SingularPotential<dim>::
return_Jacobian() const
{
    return Jac;
}



template <NematicDimension dim>
SINGULAR_POTENTIAL_EXTERNAL_LINKAGE
double SingularPotential<dim>::
return_Z() const
{
    return Z;
}



template <NematicDimension dim>
SINGULAR_POTENTIAL_EXTERNAL_LINKAGE
unsigned int SingularPotential<dim>::
invert_Q(const Eigen::Vector<double, QTensorShape<dim>::n_degrees_of_freedom> &Q_in)
{
    // TODO: add flag to reinitialize SingularPotential or not
    // TODO: figure out how to reuse Jacobian easily
    initializeInversion(Q_in);

    // Run Newton's method until residual < tolerance or reach max iterations
    for (unsigned int iter = 0; iter < maximum_iterations; ++iter)
    {
        if (Res.norm() < tolerance)
            return iter;

        this->updateVariation();
        Lambda -= damping_parameter * dLambda;
        this->updateResJac();
    }

    throw std::runtime_error("Did not calculate singular potential in under maximum iterations");
}



template<NematicDimension dim>
void SingularPotential<dim>::
initializeInversion(const Eigen::Vector<double, QTensorShape<dim>::n_degrees_of_freedom> &Q_in)
{
    inverted = false;

    m = Q_in + (1.0 / 3.0) * delta_vec;
    Lambda.setZero();
    Res = -Q_in; // can explicitly compute for Lambda = 0

    // for Jacobian, compute 2/15 on diagonal, 0 elsewhere for Lambda = 0
    for (unsigned int i = 0; i < QTensorShape<dim>::n_degrees_of_freedom; ++i)
        for (unsigned int j = 0; j < QTensorShape<dim>::n_degrees_of_freedom; ++j)
        {
            if (i == j)
                Jac(i, j) = 2.0 / 15.0;
            else
                Jac(i, j) = 0;
        }
}



template<>
void SingularPotential<NematicDimension::full_3D>::
updateResJac()
{
    Z = 0;
    Res.setZero();
    Jac.setZero();

    I2 = {};
    I4 = {};
    exp_lambda = 0;

	// Calculate each term in Lebedev quadrature for each integral, add to total
	// quadrature value until we've summed all terms
    for (std::size_t q = 0; q < w.size(); ++q)
    {
        exp_lambda = std::exp( Lambda(0) * (x[q] * x[q] - z[q] * z[q])
                               + Lambda(1) * (y[q] * y[q] - z[q] * z[q])
                               + 2 * Lambda(2) * x[q] * y[q]
                               + 2 * Lambda(3) * x[q] * z[q]
                               + 2 * Lambda(4) * y[q] * z[q] );

        Z += exp_lambda * w[q];

        I2[0] += x[q] * x[q] * exp_lambda * w[q];
        I2[1] += x[q] * y[q] * exp_lambda * w[q];
        I2[2] += y[q] * y[q] * exp_lambda * w[q];
        I2[3] += x[q] * z[q] * exp_lambda * w[q];
        I2[4] += y[q] * z[q] * exp_lambda * w[q];
        I2[5] += z[q] * z[q] * exp_lambda * w[q];

        I4[0] += x[q] * x[q] * x[q] * x[q] * exp_lambda * w[q];
        I4[1] += x[q] * x[q] * x[q] * y[q] * exp_lambda * w[q];
        I4[2] += x[q] * x[q] * y[q] * y[q] * exp_lambda * w[q];
        I4[3] += x[q] * y[q] * y[q] * y[q] * exp_lambda * w[q];
        I4[4] += y[q] * y[q] * y[q] * y[q] * exp_lambda * w[q];
        I4[5] += x[q] * x[q] * x[q] * z[q] * exp_lambda * w[q];
        I4[6] += x[q] * x[q] * y[q] * z[q] * exp_lambda * w[q];
        I4[7] += x[q] * y[q] * y[q] * z[q] * exp_lambda * w[q];
        I4[8] += y[q] * y[q] * y[q] * z[q] * exp_lambda * w[q];
        I4[9] += x[q] * x[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[10] += x[q] * y[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[11] += y[q] * y[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[12] += x[q] * z[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[13] += y[q] * z[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[14] += z[q] * z[q] * z[q] * z[q] * exp_lambda * w[q];
    }

    Res(0) = 1 / Z * I2[0] - m(0);
    Res(1) = 1 / Z * I2[2] - m(1);
    Res(2) = 1 / Z * I2[1] - m(2);
    Res(3) = 1 / Z * I2[3] - m(3);
    Res(4) = 1 / Z * I2[4] - m(4);
    
    Jac(0, 0) = 1 / Z * (I4[0] - I4[9] - 1 / Z * I2[0]*(I2[0] - I2[5]));
    Jac(0, 1) = 1 / Z * (I4[2] - I4[9] - 1 / Z * I2[0]*(I2[2] - I2[5]));
    Jac(0, 2) = 2 / Z * (I4[1] - 1 / Z * I2[0]*I2[1]);
    Jac(0, 3) = 2 / Z * (I4[5] - 1 / Z * I2[0]*I2[3]);
    Jac(0, 4) = 2 / Z * (I4[6] - 1 / Z * I2[0]*I2[4]);
    Jac(1, 0) = 1 / Z * (I4[2] - I4[11] - 1 / Z * I2[2]*(I2[0] - I2[5]));
    Jac(1, 1) = 1 / Z * (I4[4] - I4[11] - 1 / Z * I2[2]*(I2[2] - I2[5]));
    Jac(1, 2) = 2 / Z * (I4[3] - 1 / Z * I2[2]*I2[1]);
    Jac(1, 3) = 2 / Z * (I4[7] - 1 / Z * I2[2]*I2[3]);
    Jac(1, 4) = 2 / Z * (I4[8] - 1 / Z * I2[2]*I2[4]);
    Jac(2, 0) = 1 / Z * (I4[1] - I4[10] - 1 / Z * I2[1]*(I2[0] - I2[5]));
    Jac(2, 1) = 1 / Z * (I4[3] - I4[10] - 1 / Z * I2[1]*(I2[2] - I2[5]));
    Jac(2, 2) = 2 / Z * (I4[2] - 1 / Z * I2[1]*I2[1]);
    Jac(2, 3) = 2 / Z * (I4[6] - 1 / Z * I2[1]*I2[3]);
    Jac(2, 4) = 2 / Z * (I4[7] - 1 / Z * I2[1]*I2[4]);
    Jac(3, 0) = 1 / Z * (I4[5] - I4[12] - 1 / Z * I2[3]*(I2[0] - I2[5]));
    Jac(3, 1) = 1 / Z * (I4[7] - I4[12] - 1 / Z * I2[3]*(I2[2] - I2[5]));
    Jac(3, 2) = 2 / Z * (I4[6] - 1 / Z * I2[3]*I2[1]);
    Jac(3, 3) = 2 / Z * (I4[9] - 1 / Z * I2[3]*I2[3]);
    Jac(3, 4) = 2 / Z * (I4[10] - 1 / Z * I2[3]*I2[4]);
    Jac(4, 0) = 1 / Z * (I4[6] - I4[13] - 1 / Z * I2[4]*(I2[0] - I2[5]));
    Jac(4, 1) = 1 / Z * (I4[8] - I4[13] - 1 / Z * I2[4]*(I2[2] - I2[5]));
    Jac(4, 2) = 2 / Z * (I4[7] - 1 / Z * I2[4]*I2[1]);
    Jac(4, 3) = 2 / Z * (I4[10] - 1 / Z * I2[4]*I2[3]);
    Jac(4, 4) = 2 / Z * (I4[11] - 1 / Z * I2[4]*I2[4]);

    Jac_updated = true;
}



template<>
void SingularPotential<NematicDimension::quasi_2D>::
updateResJac()
{
    Z = 0;
    Res.setZero();
    Jac.setZero();

    I2 = {};
    I4 = {};
    exp_lambda = 0;

	// Calculate each term in Lebedev quadrature for each integral, add to total
	// quadrature value until we've summed all terms
    for (std::size_t q = 0; q < w.size(); ++q)
    {
        exp_lambda = std::exp( Lambda(0) * (x[q] * x[q] - z[q] * z[q])
                               + Lambda(1) * (y[q] * y[q] - z[q] * z[q])
                               + 2 * Lambda(2) * x[q] * y[q] );

        Z += exp_lambda * w[q];

        I2[0] += x[q] * x[q] * exp_lambda * w[q];
        I2[1] += x[q] * y[q] * exp_lambda * w[q];
        I2[2] += y[q] * y[q] * exp_lambda * w[q];
        I2[3] += z[q] * z[q] * exp_lambda * w[q];

        I4[0] += x[q] * x[q] * x[q] * x[q] * exp_lambda * w[q];
        I4[1] += x[q] * x[q] * x[q] * y[q] * exp_lambda * w[q];
        I4[2] += x[q] * x[q] * y[q] * y[q] * exp_lambda * w[q];
        I4[3] += x[q] * y[q] * y[q] * y[q] * exp_lambda * w[q];
        I4[4] += y[q] * y[q] * y[q] * y[q] * exp_lambda * w[q];
        I4[5] += x[q] * x[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[6] += x[q] * y[q] * z[q] * z[q] * exp_lambda * w[q];
        I4[7] += y[q] * y[q] * z[q] * z[q] * exp_lambda * w[q];
    }

    Res(0) = 1 / Z * I2[0] - m(0);
    Res(1) = 1 / Z * I2[2] - m(1);
    Res(2) = 1 / Z * I2[1] - m(2); 

    Jac(0, 0) = 1 / Z * (I4[0] - I4[5] - 1 / Z * I2[0]*(I2[0] - I2[3]));
    Jac(0, 1) = 1 / Z * (I4[2] - I4[5] - 1 / Z * I2[0]*(I2[2] - I2[3]));
    Jac(0, 2) = 2 / Z * (I4[1] - 1 / Z * I2[0]*I2[1]);
    Jac(1, 0) = 1 / Z * (I4[2] - I4[7] - 1 / Z * I2[2]*(I2[0] - I2[3]));
    Jac(1, 1) = 1 / Z * (I4[4] - I4[7] - 1 / Z * I2[2]*(I2[2] - I2[3]));
    Jac(1, 2) = 2 / Z * (I4[3] - 1 / Z * I2[2]*I2[1]);
    Jac(2, 0) = 1 / Z * (I4[1] - I4[6] - 1 / Z * I2[1]*(I2[0] - I2[3]));
    Jac(2, 1) = 1 / Z * (I4[3] - I4[6] - 1 / Z * I2[1]*(I2[2] - I2[3]));
    Jac(2, 2) = 2 / Z * (I4[2] - 1 / Z * I2[1]*I2[1]);

    Jac_updated = true;
}



template <NematicDimension dim>
void SingularPotential<dim>::
updateVariation()
{
    dLambda = Jac.householderQr().solve(Res);
    Jac_updated = false; // Can't use Jac when it's factorized
}

// template class SingularPotential<NematicDimension::full_2D>;
template class SingularPotential<NematicDimension::quasi_2D>;
template class SingularPotential<NematicDimension::full_3D>;

} // ball_majumdar_singular_potential
