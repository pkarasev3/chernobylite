project(area_tracker_cpp)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_VERBOSE_MAKEFILE CACHE bool true)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)

set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE PATH
  "Single Directory for all Executables")
set( LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE PATH
        "Single Directory for all Libraries")

set(CMAKE_CXX_FLAGS "-Wno-long-long -Wno-variadic-macros -Wall -pedantic -pipe -msse4.2" CACHE STRING "compiler flags always used" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g3" CACHE STRING "extra compiler flags for debug")
set(CMAKE_C_FLAGS "-Wno-long-long -std=c99 -Wall -pedantic -pipe -msse4.2" CACHE STRING "compiler flags always used" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g3" CACHE STRING "extra compiler flags for debug" )
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "extra compiler flags for release")
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "configuration (default: Debug)" )

find_package(OpenCV  REQUIRED)

add_executable( match_images_opencv  match_images_opencv.cpp )
target_link_libraries( match_images_opencv ${OpenCV_LIBS} )

add_executable( match_images_opencv_display  match_images_opencv_display.cpp )
target_link_libraries( match_images_opencv_display ${OpenCV_LIBS} )


message(STATUS "    Compiler:                  ${CMAKE_COMPILER}")
message(STATUS "    C++ flags (Release):       ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "    C++ flags (Debug):         ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "OpenCV Version is ${OpenCV_VERSION} make sure you are using latest from trunk")


