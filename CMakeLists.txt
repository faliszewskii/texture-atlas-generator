cmake_minimum_required(VERSION 3.28)
project(texture_atlas_generator)
set(LIBRARY_NAME texture_atlas_generator)

set(CMAKE_CXX_STANDARD 26)

# Texture atlas generator
set(TAG_SRC
        src/model/Triangulation.h
        src/model/Mesh.h
        src/alg/LSCM.cpp
        src/alg/LSCM.h
        src/model/UvMap.h
)

add_library(${LIBRARY_NAME} STATIC ${TAG_SRC})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (${LIBRARY_NAME} Eigen3::Eigen)
target_link_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/dep/eigen)
target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/dep/eigen)
