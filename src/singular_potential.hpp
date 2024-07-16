#ifndef SINGULAR_POTENTIAL_HPP
#define SINGULAR_POTENTIAL_HPP

#include "preprocessor.hpp"
#include "Q_tensor_shape.hpp"

#include <cmath>

#include "../eigen/Eigen/Dense"
#include "../lebedev-quadrature/src/lebedev_quadrature.hpp"

namespace ball_majumdar_singular_potential
{

template <NematicDimension dim>
class SingularPotential
{
public:
    using mat = Eigen::Matrix<double, 
                              QTensorShape<dim>::n_degrees_of_freedom,
                              QTensorShape<dim>::n_degrees_of_freedom>;
    using vec = Eigen::Vector<double, QTensorShape<dim>::n_degrees_of_freedom>;
    using point = Eigen::Vector<double, QTensorShape<dim>::matrix_dimension>;

    SingularPotential(const unsigned int lebedev_order,
                      const double damping_parameter,
                      const double tolerance,
                      const unsigned int maximum_iterations);

    unsigned int invert_Q(const vec &Q_in);
    double return_Z() const;
    vec return_Lambda() const;
    mat return_Jacobian() const;

private:
    void initializeInversion(const vec &Q_in);
    void updateResJac();
    void updateVariation();

    const double damping_parameter;
    const double tolerance;
    const unsigned int maximum_iterations;

    bool inverted = false;
    bool Jac_updated = false;

    vec delta_vec;
    vec m;
    vec Q;
    vec Lambda;

    vec Res;
    mat Jac;
    vec dLambda;
    double Z = 0;

    double exp_lambda = 0;
    std::array<double, QTensorShape<dim>::degree_2_integrals> I2 = {0};
    std::array<double, QTensorShape<dim>::degree_4_integrals> I4 = {0};

    lebedev::QuadraturePoints quadrature_points;

    const std::vector<double> &x;
    const std::vector<double> &y;
    const std::vector<double> &z;
    const std::vector<double> &w;
};

} // ball_majumdar_singular_potential

#endif
