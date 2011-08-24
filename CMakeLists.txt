cmake_minimum_required(VERSION 2.8)
project(chernobylite)

set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "The build type.")

find_package( Boost COMPONENTS filesystem program_options thread )
find_package( OpenCV REQUIRED )

link_directories( ${Boost_LIBRARY_DIRS} )
include_directories( ${Boost_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} )
message(STATUS "boost include dirs:  " ${Boost_INCLUDE_DIRS} )
message(STATUS "boost library dirs:  " ${Boost_LIBRARY_DIRS} )

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

add_definitions("-Werror -fno-strict-aliasing -Wall -Wl,--no-undefined")
include_directories(${chernobylite_SOURCE_DIR}/include)

#put impl of all the helpers into library.
add_library(chernobylite SHARED src/chernobylite.cpp)
target_link_libraries(chernobylite ${OpenCV_LIBS} ${Boost_LIBRARIES})
# add libraries and such.

add_subdirectory(unittest)
