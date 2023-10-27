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


set(SOLARSHADING_TARGET_NAME feelpp_solarshading)
add_library(${SOLARSHADING_TARGET_NAME} INTERFACE)
target_include_directories(${SOLARSHADING_TARGET_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )
target_link_libraries(${SOLARSHADING_TARGET_NAME} INTERFACE Feelpp::feelpp)

feelpp_add_application(example_shading_masks_surface_mesh SRCS example_shading_masks.cpp DEFS FEELPP_TOP_DIM=2 LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME} TESTS NO_MPI_TEST INSTALL)
feelpp_add_application(example_shading_masks_volume_mesh SRCS example_shading_masks.cpp DEFS FEELPP_TOP_DIM=3 LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME} TESTS NO_MPI_TEST INSTALL)

feelpp_add_application(shading_mask_ai_dataset SRCS shading_mask_ai_dataset.cpp LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME} INSTALL)


feelpp_add_testcase(cases)
