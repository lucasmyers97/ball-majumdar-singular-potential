#ifndef Q_TENSOR_SHAPE_HPP
#define Q_TENSOR_SHAPE_HPP

#include "preprocessor.hpp"

#include <array>

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

    using vec = std::array<int, n_degrees_of_freedom>;
    static constexpr vec dof_tensor_row_idx = { {0, 0} };
    static constexpr vec dof_tensor_col_idx = { {0, 1} };
    static constexpr vec delta_vec = { {1, 0} };
};



template <>
struct QTensorShape<NematicDimension::quasi_2D>
{
    static constexpr unsigned int matrix_dimension = 3;
    static constexpr unsigned int n_degrees_of_freedom = 3;

    using vec = std::array<int, n_degrees_of_freedom>;
    static constexpr vec dof_tensor_row_idx = { {0, 1, 0} };
    static constexpr vec dof_tensor_col_idx = { {0, 1, 1} };
    static constexpr vec delta_vec = { {1, 1, 0} };
};



template <>
struct QTensorShape<NematicDimension::full_3D>
{
    static constexpr unsigned int matrix_dimension = 3;
    static constexpr unsigned int n_degrees_of_freedom = 5;

    using vec = std::array<int, n_degrees_of_freedom>;
    static constexpr vec dof_tensor_row_idx = { {0, 1, 0, 0, 1} };
    static constexpr vec dof_tensor_col_idx = { {0, 1, 1, 2, 2} };
    static constexpr vec delta_vec = { {1, 1, 0, 0, 0} };
};

} // namespace ball_majumdar_singular_potential

#endif
