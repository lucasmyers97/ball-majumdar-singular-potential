add_executable(matrix_inversion matrix_inversion.cpp)
target_link_libraries(matrix_inversion
    PRIVATE
    Eigen3::Eigen)

install(TARGETS
    matrix_inversion
    DESTINATION bin/examples)
