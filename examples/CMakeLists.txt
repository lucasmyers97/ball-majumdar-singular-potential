add_executable(matrix_inversion matrix_inversion.cpp)
target_link_libraries(matrix_inversion
    PRIVATE
    Eigen)

add_executable(print_Q_tensor_shape print_Q_tensor_shape.cpp)
target_link_libraries(print_Q_tensor_shape
    PRIVATE
    Q_tensor_shape)

add_executable(lebedev_quadrature_test lebedev_quadrature_test.cpp)
target_link_libraries(lebedev_quadrature_test
    PRIVATE
    lebedev_quadrature)

add_executable(singular_potential_test singular_potential_test.cpp)
target_link_libraries(singular_potential_test
    PRIVATE
    singular_potential
    Q_tensor_shape)

install(TARGETS
    matrix_inversion
    print_Q_tensor_shape
    lebedev_quadrature_test
    singular_potential_test
    DESTINATION bin/examples)
