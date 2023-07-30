add_executable(matrix_inversion matrix_inversion.cpp)
target_link_libraries(matrix_inversion
    PRIVATE
    Eigen)

add_executable(print_Q_tensor_shape print_Q_tensor_shape.cpp)
target_link_libraries(print_Q_tensor_shape
    PRIVATE
    Q_tensor_shape)

install(TARGETS
    matrix_inversion
    print_Q_tensor_shape
    DESTINATION bin/examples)
