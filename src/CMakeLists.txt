add_library(Q_tensor_shape INTERFACE)
target_include_directories(Q_tensor_shape INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(singular_potential INTERFACE)
target_include_directories(singular_potential
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR})
# target_link_libraries(singular_potential
#     INTERFACE
#     Q_tensor_shape
#     Eigen
#     lebedev_quadrature)
target_link_libraries(singular_potential
    INTERFACE
    Q_tensor_shape)
