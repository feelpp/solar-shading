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


feelpp_add_application(example_ShadingMasks SRCS example_ShadingMasks.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


feelpp_add_testcase(cases)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BELOW THE IN-APP BENCHMARKING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Add the benchmarkPCG target
add_executable(example_ShadingMasks_comparison example_ShadingMasks_comparison.cpp)

feelpp_add_application(example_ShadingMasks_comparison SRCS example_ShadingMasks_comparison.cpp LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)


target_include_directories(example_ShadingMasks_comparison PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigen/Eigen
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
)

target_link_libraries(example_ShadingMasks_comparison PRIVATE
    Feelpp::feelpp
)

# Add compilation flags
target_compile_options(example_ShadingMasks_comparison PRIVATE
    -O3
    -march=native
    -funroll-loops
)
