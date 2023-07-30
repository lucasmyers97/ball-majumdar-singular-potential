#ifndef Q_TENSOR_SHAPE_HPP
#define Q_TENSOR_SHAPE_HPP

namespace ball_majumdar_singular_potential
{

enum class NematicDimension
{
    full_2D,
    quasi_2D,
    full_3D,
};



template <NematicDimension dim>
struct QTensorShape
{};



template <>
struct QTensorShape<NematicDimension::full_2D>
{
    static constexpr unsigned int matrix_dimension = 2;
    static constexpr unsigned int n_degrees_of_freedom = 2;
};



template <>
struct QTensorShape<NematicDimension::quasi_2D>
{
    static constexpr unsigned int matrix_dimension = 3;
    static constexpr unsigned int n_degrees_of_freedom = 3;
};



template <>
struct QTensorShape<NematicDimension::full_3D>
{
    static constexpr unsigned int matrix_dimension = 3;
    static constexpr unsigned int n_degrees_of_freedom = 5;
};

} // namespace ball_majumdar_singular_potential

#endif
