#include <lebedev_quadrature.hpp>

#include <iostream>
#include <cmath>

double xyz_squared(double x, double y, double z)
{
    return x*x * y*y * z*z;
}



int main()
{
    auto quad_order = lebedev::QuadratureOrder::order_590;
    auto quad_points = lebedev::QuadraturePoints(quad_order);

    double quadrature_value = quad_points.evaluate_spherical_integral(xyz_squared);

    std::cout << "Integral is: " << 4 * M_PI * quadrature_value << "\n";
}
