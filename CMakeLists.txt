cmake_minimum_required(VERSION 3.16)

project(nbody_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(MPI REQUIRED)
include(CTest)

add_library(nbody_core
    src/simulation.cpp
    src/particle_io.cpp
    src/distribution.cpp
    src/octree.cpp
)

target_include_directories(nbody_core PUBLIC src)

add_executable(nbody src/main.cpp)

target_link_libraries(nbody PRIVATE nbody_core MPI::MPI_CXX)

target_compile_options(nbody_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)
target_compile_options(nbody PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)

if(BUILD_TESTING)
    add_executable(nbody_tests tests/test_core.cpp)
    target_link_libraries(nbody_tests PRIVATE nbody_core)
    target_compile_options(nbody_tests PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
    )
    add_test(NAME nbody_core COMMAND nbody_tests)
endif()
