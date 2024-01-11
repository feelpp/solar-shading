###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 07 Apr 2019
#
#  Copyright (C) 2019 Feel++ Consortium
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



#BEGIN::SPECX
# will need to be replaced
file(GLOB_RECURSE Specx_list 
    "/home/u2/lemoine/SpecxProjects/specxpp/*.hpp"
)
if (TRUE)
    foreach(file_path ${Specx_list})
        message(${file_path})
    endforeach()
endif()

MACRO(HEADER_DIRECTORIES return_list)
    FILE(GLOB_RECURSE new_list /home/u2/lemoine/SpecxProjects/specxpp/Src/*.hpp)
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

HEADER_DIRECTORIES(header_specx_dir_list)

list(LENGTH header_specx_dir_list header_dir_list_count)
message(STATUS "[INFO] Found ${header_dir_list_count} header directories.")

if (TRUE)
    foreach(dir_path ${header_specx_dir_list})
        message(${dir_path})
    endforeach()
endif()

include_directories(/home/u2/lemoine/SpecxProjects/specxpp/build/Src) 
include_directories(/home/u2/lemoine/SpecxProjects/specxpp/Src) 
#END::SPECX



set(SOLARSHADING_TARGET_NAME_SPECX feelpp_solarshading_specx)
add_library(${SOLARSHADING_TARGET_NAME_SPECX} INTERFACE)


# will need to be replaced
target_include_directories(${SOLARSHADING_TARGET_NAME_SPECX} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
    ${header_specx_dir_list}
    #Feelpp::feelpp_specx
  )

target_link_libraries(${SOLARSHADING_TARGET_NAME_SPECX} 
    INTERFACE 
        Feelpp::feelpp
        Feelpp::feelpp_specx
)




if ( FEELPP_HAS_SPECX )
    MESSAGE(***********+++++++++ S P E C X in FEELPP ++++++++++++*******************)
endif()

if(TRUE)
    #feelpp_add_application( SM_Specx_3 SRCS Specx_Solar_Shading.cpp DEFS FEELPP_TOP_DIM=3 LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME_SPECX} Feelpp::feelpp_specx )
    feelpp_add_application(example_shading_masks_volume_mesh SRCS example_shading_masks.cpp DEFS FEELPP_TOP_DIM=3 LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME_SPECX} Feelpp::feelpp_specx)
endif()

if(TRUE)
    #feelpp_add_application( SM_Specx_2 SRCS Specx_Solar_Shading.cpp DEFS FEELPP_TOP_DIM=2 LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME_SPECX} Feelpp::feelpp_specx )    
    feelpp_add_application(example_shading_masks_surface_mesh SRCS example_shading_masks.cpp DEFS FEELPP_TOP_DIM=2 LINK_LIBRARIES ${SOLARSHADING_TARGET_NAME_SPECX} Feelpp::feelpp_specx)
endif()

if(TRUE)
    feelpp_add_testcase(cases)
endif()


