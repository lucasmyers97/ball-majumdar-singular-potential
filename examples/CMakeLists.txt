add_executable(matrix_inversion matrix_inversion.cpp)

add_executable(print_Q_tensor_shape print_Q_tensor_shape.cpp)
target_link_libraries(print_Q_tensor_shape
    PRIVATE
    Q_tensor_shape)

add_executable(lebedev_quadrature_test lebedev_quadrature_test.cpp)

add_executable(singular_potential_test singular_potential_test.cpp)
target_link_libraries(singular_potential_test
    PRIVATE
    singular_potential
    Q_tensor_shape)

add_library(singular_potential_implementation
    singular_potential_implementation.cpp)
target_link_libraries(singular_potential_implementation
    PUBLIC
    singular_potential)
target_include_directories(singular_potential_implementation
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(singular_potential_implementation_test singular_potential_implementation_test.cpp)
target_link_libraries(singular_potential_implementation_test
    PRIVATE
    singular_potential_implementation)

install(TARGETS
    matrix_inversion
    print_Q_tensor_shape
    lebedev_quadrature_test
    singular_potential_test
    singular_potential_implementation_test
    DESTINATION bin/examples)
