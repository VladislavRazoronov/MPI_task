cmake_minimum_required(VERSION 3.17)
project(MPI_project)

set(CMAKE_CXX_STANDARD 14)

FIND_PACKAGE(Boost 1.4 REQUIRED)
find_package(ImageMagick COMPONENTS Magick++ convert REQUIRED)
include_directories(MPI_project ${Boost_INCLUDE_DIRS})
include_directories(MPI_project ${ImageMagick_INCLUDE_DIRS})

add_executable(MPI_project main.cpp Functions.cpp Functions.h)
add_definitions( -DMAGICKCORE_QUANTUM_DEPTH=16 )
add_definitions( -DMAGICKCORE_HDRI_ENABLE=0 )
target_link_libraries(MPI_project ${Boost_LIBRARIES})
target_link_libraries(MPI_project ${ImageMagick_LIBRARIES})
