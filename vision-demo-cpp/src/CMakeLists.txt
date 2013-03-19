cmake_minimum_required(VERSION 2.8)
project(VisionDemoChernobylite)

enable_testing()

set(CMAKE_VERBOSE_MAKEFILE CACHE bool true)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)

set( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE PATH
  "Single Directory for all Executables")
set( LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE PATH
        "Single Directory for all Libraries")

set(CMAKE_CXX_FLAGS "-std=c++0x -Wno-long-long -Wno-variadic-macros -Wall -pedantic -pipe -msse4.2" CACHE STRING "compiler flags always used" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g3" CACHE STRING "extra compiler flags for debug")
set(CMAKE_C_FLAGS "-Wno-long-long  -Wall -pedantic -pipe -msse4.2" CACHE STRING "compiler flags always used" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g3" CACHE STRING "extra compiler flags for debug" )
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "extra compiler flags for release")
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "configuration (default: Debug)" )

find_package(OpenCV  REQUIRED)
find_package( Boost COMPONENTS filesystem program_options )
find_package( OpenMP )


link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )

add_executable( match_images_opencv  apps/match_images_opencv.cpp )
target_link_libraries( match_images_opencv ${OpenCV_LIBS} )

add_executable( match_images_opencv_display  apps/match_images_opencv_display.cpp )
target_link_libraries( match_images_opencv_display ${OpenCV_LIBS} )


message(STATUS "    Compiler:                  ${CMAKE_COMPILER}")
message(STATUS "    C++ flags (Release):       ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "    C++ flags (Debug):         ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "OpenCV Version is ${OpenCV_VERSION} make sure you are using a sane version")


#add_executable( WebcamViewer  apps/WebcamViewer.cpp )
#target_link_libraries( WebcamViewer ${OpenCV_LIBS} ${Boost_LIBRARIES} )


##################################################### OpenGL stuff
find_package(OpenGL REQUIRED)
find_package(GLUT)
set(CORE3DLIBS X11 ${GLUT_LIBRARY} ${OPENGL_LIBRARY} GLEW GLU m)

file(GLOB glToolSrc GLTools/src/*.cpp ) #GLTools/src/*.c
file(GLOB glToolHdr GLTools/include/*.h)
message(STATUS "   ${OPENGL_LIBRARY},  glTools src & hdr:         ${glToolSrc} ${glToolHdr}")

include_directories(${CMAKE_SOURCE_DIR}/GLTools/include)
add_executable(GSExplode  apps/GSExplode.cpp  ${glToolSrc} ${glToolHdr})
target_link_libraries(GSExplode  ${CORE3DLIBS} ${OpenCV_LIBS} ${Boost_LIBRARIES})
