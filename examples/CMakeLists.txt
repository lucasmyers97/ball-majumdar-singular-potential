add_executable(matrix_inversion matrix_inversion.cpp)
target_link_libraries(matrix_inversion
    PRIVATE
    Eigen)

install(TARGETS
    matrix_inversion
    DESTINATION bin/examples)
