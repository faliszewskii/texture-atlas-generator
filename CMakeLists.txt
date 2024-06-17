cmake_minimum_required(VERSION 3.28)
project(texture_atlas_generator)

set(CMAKE_CXX_STANDARD 26)

add_library(texture_atlas_generator STATIC library.cpp)

include_directories(${CMAKE_SOURCE_DIR}/dep/eigen)