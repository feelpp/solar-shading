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


feelpp_add_application(example_ShadingMasks SRCS example_ShadingMasks.cpp spa.c LINK_LIBRARIES Feelpp::feelpp TESTS NO_MPI_TEST INSTALL)

feelpp_add_testcase(cases)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BELOW THE IN-APP BENCHMARKING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Benchmark with compilation flags -Ofast -march=native

add_executable(example_ShadingMasks_comparison1 example_ShadingMasks_comparison.cpp)
target_include_directories(example_ShadingMasks_comparison1 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)
target_link_libraries(example_ShadingMasks_comparison1 PRIVATE
    Feelpp::feelpp
)
target_compile_options(example_ShadingMasks_comparison1 PRIVATE
    -Ofast
    -Wno-date-time
)

# Benchmark with compilation flags -03 -march=native

add_executable(example_ShadingMasks_comparison2 example_ShadingMasks_comparison.cpp)
target_include_directories(example_ShadingMasks_comparison2 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)
target_link_libraries(example_ShadingMasks_comparison2 PRIVATE
    Feelpp::feelpp
)
target_compile_options(example_ShadingMasks_comparison2 PRIVATE
    -O3
    -Wno-date-time
)

# Benchmark with compilation flags -O2 -march=native

add_executable(example_ShadingMasks_comparison3 example_ShadingMasks_comparison.cpp)

target_include_directories(example_ShadingMasks_comparison3 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)
target_link_libraries(example_ShadingMasks_comparison3 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison3 PRIVATE
    -O2
    -Wno-date-time
)

# Benchmark with compilation flags -Ofast -march=native -funroll-loops

add_executable(example_ShadingMasks_comparison4 example_ShadingMasks_comparison.cpp)

target_include_directories(example_ShadingMasks_comparison4 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)

target_link_libraries(example_ShadingMasks_comparison4 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison4 PRIVATE
    -Ofast
    -funroll-loops
    -Wno-date-time
)

# Benchmark with compilation flags -O3 -march=native -funroll-loops

add_executable(example_ShadingMasks_comparison5 example_ShadingMasks_comparison.cpp)

target_include_directories(example_ShadingMasks_comparison5 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)

target_link_libraries(example_ShadingMasks_comparison5 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison5 PRIVATE
    -O3
    -funroll-loops
    -Wno-date-time
)

# Benchmark with compilation flags -O2 -march=native -funroll-loops

add_executable(example_ShadingMasks_comparison6 example_ShadingMasks_comparison.cpp)
target_include_directories(example_ShadingMasks_comparison6 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)

target_link_libraries(example_ShadingMasks_comparison6 PUBLIC
    Feelpp::feelpp
)

install(TARGETS example_ShadingMasks_comparison6 RUNTIME DESTINATION bin COMPONENT Bin)

target_compile_options(example_ShadingMasks_comparison6 PRIVATE
    -O2
    -funroll-loops
    -Wno-date-time
)

# Benchmark with compilation flags -Ofast -march=native -funroll-loops -funsafe-math-optimizations

add_executable(example_ShadingMasks_comparison7 example_ShadingMasks_comparison.cpp)

target_include_directories(example_ShadingMasks_comparison7 PRIVATE
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/pcg-cpp/include
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/Xoshiro-cpp
    ${CMAKE_SOURCE_DIR}/benchmark/extlibs/eigenrand
)

target_link_libraries(example_ShadingMasks_comparison7 PRIVATE
    Feelpp::feelpp
)

target_compile_options(example_ShadingMasks_comparison7 PRIVATE
    -Ofast
    -funroll-loops
    -funsafe-math-optimizations
    -Wno-date-time
)