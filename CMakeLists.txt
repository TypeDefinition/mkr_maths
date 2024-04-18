cmake_minimum_required(VERSION 3.23.2)
project(mkr_maths_lib)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-std=c++20 -O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Source Files
set(SRC_DIR "src")
file(GLOB_RECURSE SRC_FILES LIST_DIRECTORIES true CONFIGURE_DEPENDS
        "${SRC_DIR}/*.h"
        "${SRC_DIR}/*.c"
        "${SRC_DIR}/*.hpp"
        "${SRC_DIR}/*.cpp")
add_library(${PROJECT_NAME} ${SRC_FILES})

# Target
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_DIR})

# Test
enable_testing()
add_subdirectory(test)