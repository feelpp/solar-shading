###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2020-06-03
#
#  Copyright (C) 2020 Cemosis
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

cmake_minimum_required(VERSION 3.21.0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USING THE GPU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project(ShadingMasksGPU LANGUAGES CXX CUDA)

set(CMAKE DISABLE_IN_SOURCE_BUILD ON)

# If you want to use clang as your C++ compiler
set(CMAKE_CXX_COMPILER clang++)

find_package(Feel++ REQUIRED)

# Add executable target with source files
add_executable(ShadingMasksGPU example_ShadingMasks_GPU.cpp Bvh_GPU.cu Bvh_GPU.cuh)

# Specify the C++ standard
set_target_properties(ShadingMasksGPU PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_EXTENSIONS OFF CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES "86")

# Add the Feel++ libraries to the link libraries
target_link_libraries(ShadingMasksGPU PUBLIC Feelpp::feelpp)

target_compile_options(ShadingMasksGPU PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_86,code=sm_86>)

target_include_directories(ShadingMasksGPU PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}/src)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ NORMAL IMPLEMENTATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

feelpp_add_application(example_ShadingMasks SRCS example_ShadingMasks.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)

feelpp_add_testcase(cases)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BELOW THE IN-APP BENCHMARKING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Benchmark with compilation flags -Ofast -march=native

add_executable(example_ShadingMasks_comparison1 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison1 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison1 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison1 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison1 PRIVATE
    -Ofast
    -march=native
    -Wno-date-time
)

# Benchmark with compilation flags -03 -march=native

add_executable(example_ShadingMasks_comparison2 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison2 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison2 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison2 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison2 PRIVATE
    -O3
    -march=native
    -Wno-date-time
)

# Benchmark with compilation flags -O2 -march=native

add_executable(example_ShadingMasks_comparison3 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison3 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison3 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison3 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison3 PRIVATE
    -O2
    -march=native
    -Wno-date-time
)

# Benchmark with compilation flags -Ofast -march=native -funroll-loops

add_executable(example_ShadingMasks_comparison4 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison4 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison4 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison4 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison4 PRIVATE
    -Ofast
    -march=native
    -funroll-loops
    -Wno-date-time
)

# Benchmark with compilation flags -O3 -march=native -funroll-loops

add_executable(example_ShadingMasks_comparison5 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison5 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison5 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison5 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison5 PRIVATE
    -O3
    -march=native
    -funroll-loops
    -Wno-date-time
)

# Benchmark with compilation flags -O2 -march=native -funroll-loops

add_executable(example_ShadingMasks_comparison6 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison6 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison6 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison6 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison6 PRIVATE
    -O2
    -march=native
    -funroll-loops
    -Wno-date-time
)

# Benchmark with compilation flags -Ofast -march=native -funroll-loops -funsafe-math-optimizations

add_executable(example_ShadingMasks_comparison7 example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison7 SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison7 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison7 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison7 PRIVATE
    -Ofast
    -march=native
    -funroll-loops
    -funsafe-math-optimizations
    -Wno-date-time
)