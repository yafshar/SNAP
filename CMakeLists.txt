#
# CMakeLists.txt
#
# LGPL Version 2.1 HEADER START
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
#
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301  USA
#
# LGPL Version 2.1 HEADER END
#

#
# Copyright (c) 2013--2020, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Richard Berger
#    Christoph Junghans
#    Ryan S. Elliott
#    Yaser Afshar
#


# Required preamble

cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API 2.0 REQUIRED CONFIG)
if(NOT TARGET kim-api)
  enable_testing()
  project("${KIM_API_PROJECT_NAME}" VERSION "${KIM_API_VERSION}"
    LANGUAGES CXX C Fortran)
endif()

# End preamble

# This model driver requires C++11
set(CMAKE_CXX_STANDARD 11)

set(MODEL_DRIVER_NAME "SNAP__MD_536750310735_001")

add_kim_api_model_driver_library(
  NAME                    ${MODEL_DRIVER_NAME}
  CREATE_ROUTINE_NAME     "model_driver_create"
  CREATE_ROUTINE_LANGUAGE "cpp"
  )

target_sources(${MODEL_DRIVER_NAME} PRIVATE
  SNAP.cpp
  SNAPImplementation.cpp
  SNA.cpp
  ZBL.cpp
  TABLE.cpp
  )
