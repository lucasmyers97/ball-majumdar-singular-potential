#include "Q_tensor_shape.hpp"

#include <iostream>

using namespace ball_majumdar_singular_potential;

int main()
{
    using QShape2D = QTensorShape<NematicDimension::full_2D>;
    using QShapeQuasi2D = QTensorShape<NematicDimension::quasi_2D>;
    using QShape3D = QTensorShape<NematicDimension::full_3D>;

    std::cout << "2D matrix dimension: "
              << QShape2D::matrix_dimension
              << "\n" << "2D degrees of freedom: "
              << QShape2D::n_degrees_of_freedom << "\n";

    std::cout << "quasi-2D matrix dimension: "
              << QShapeQuasi2D::matrix_dimension
              << "\n" << "quasi-2D degrees of freedom: "
              << QShapeQuasi2D::n_degrees_of_freedom << "\n";

    std::cout << "3D matrix dimension: "
              << QShape3D::matrix_dimension
              << "\n" << "3D degrees of freedom: "
              << QShape3D::n_degrees_of_freedom << "\n";

    return 0;
}
