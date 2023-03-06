# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

project(ball-majumdar-singular-potential VERSION 0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_COMPILER "clang++-13")
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "Compiler")
add_compile_options(-Wall -Wextra -pedantic)
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# library source code
add_subdirectory(src)
