# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

project(ball-majumdar-singular-potential VERSION 0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_COMPILER "clang++-13")
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "Compiler")
add_compile_options(-Wall -Wextra -pedantic)
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# have eigen installed as submodule, header only
# add_library(Eigen INTERFACE)
# target_include_directories(Eigen INTERFACE ${CMAKE_SOURCE_DIR}/eigen)
# 
# # have lebedev-quadrature installed as submodule, header only
# add_library(lebedev_quadrature INTERFACE)
# target_include_directories(lebedev_quadrature INTERFACE ${CMAKE_SOURCE_DIR}/lebedev-quadrature/src)

# library source code
add_subdirectory(src)

# examples for how to use the library
add_subdirectory(examples)
