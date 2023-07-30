#ifndef SINGULAR_POTENTIAL_HPP
#define SINGULAR_POTENTIAL_HPP

#include "Q_tensor_shape.hpp"

#include <Eigen/Dense>
#include <lebedev_quadrature.hpp>

namespace ball_majumdar_singular_potential
{

template <NematicDimension dim>
class SingularPotential
{
private:
    using mat = Eigen::Matrix<double, 
                              QTensorShape<dim>::n_degrees_of_freedom,
                              QTensorShape<dim>::n_degrees_of_freedom>;
    using vec = Eigen::Vector<double, QTensorShape<dim>::n_degrees_of_freedom>;
    using point = Eigen::Vector<double, QTensorShape<dim>::matrix_dimension>;

public:

    SingularPotential(const unsigned int lebedev_order,
                      const double damping_parameter,
                      const double tolerance,
                      const unsigned int maximum_iterations);

    unsigned int invert_Q(const vec &Q_in);
    double returnZ() const;
    void returnLambda(vec &outLambda) const;
    void returnJac(mat &outJac);

private:
    void initializeInversion(const vec &Q_in);
    void updateResJac();
    void updateVariation();

    double lambdaSum(const point &x) const;
    double calcInt1Term(const double exp_lambda,
                        const int quad_idx, const int i_m, const int j_m) const;
    double calcInt2Term(const double exp_lambda,
                        const int quad_idx, const int i_m, const int j_m,
                        const int i_n, const int j_n) const;
    double calcInt3Term(const double exp_lambda,
                        const int quad_idx,
                        const int i_m, const int j_m,
                        const int i_n, const int j_n) const;
    double calcInt4Term(const double exp_lambda,
                        const int quad_idx, const int i_m, const int j_m) const;

    const double damping_parameter;
    const double tolerance;
    const unsigned int maximum_iterations;

    bool inverted = false;
    bool Jac_updated = false;

    vec Q;
    vec Lambda;

    vec Res;
    mat Jac;
    vec dLambda;
    double Z = 0;

    using int_vec = std::array<double, QTensorShape<dim>::n_degrees_of_freedom>;
    using int_mat = std::array<int_vec, QTensorShape<dim>::n_degrees_of_freedom>;

    int_vec int1 = {0};
    int_mat int2 = {{0}};
    int_mat int3 = {{0}};
    int_vec int4 = {0};

    lebedev::QuadraturePoints quadrature_points;
};

} // ball_majumdar_singular_potential

#endif
