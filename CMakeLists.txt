# CMake Version
cmake_minimum_required(VERSION 3.23.2)
# Project Name
project(mkr_maths_lib)

# Includes
include(FetchContent)

# C++ Version
set(CMAKE_CXX_STANDARD 20)

# Compiler Flags
set(CMAKE_CXX_FLAGS "-std=c++20 -O3")

# mkr Engine
set(SRC_DIR "src/")
file(GLOB_RECURSE SRC_FILES LIST_DIRECTORIES true CONFIGURE_DEPENDS
        "${SRC_DIR}/*.h"
        "${SRC_DIR}/*.c"
        "${SRC_DIR}/*.hpp"
        "${SRC_DIR}/*.cpp")
add_library(${PROJECT_NAME} ${SRC_FILES})

# External Dependencies
FetchContent_Declare(mkr_common_lib GIT_REPOSITORY https://github.com/TypeDefinition/mkr_common_lib.git GIT_TAG main)
FetchContent_MakeAvailable(mkr_common_lib)

# Target
target_include_directories(${PROJECT_NAME} PUBLIC ${SRC_DIR})
target_link_libraries(${PROJECT_NAME} mkr_common_lib)