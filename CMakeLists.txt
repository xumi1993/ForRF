cmake_minimum_required(VERSION 3.10)
project(ForRF LANGUAGES Fortran CUDA CXX)


# Set Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

enable_testing()
add_subdirectory(src)

