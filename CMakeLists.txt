cmake_minimum_required(VERSION 3.18)
project(towed_sonar_eq
    VERSION 1.0.0
    DESCRIPTION "核潜艇拖曳线列阵声呐时域基带信号处理与水声通信还原内核"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")

option(BUILD_PYTHON_BINDINGS "Build Python bindings with pybind11" ON)
option(BUILD_TESTS "Build unit tests" ON)

add_library(sonar_core STATIC
    src/hydrophone_unpacker.cpp
    src/cma_equalizer.cpp
    src/acoustic_ray_tracer.cpp
    src/doppler_resampler.cpp
    src/spectrum_engine.cpp
)

target_include_directories(sonar_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(sonar_core PUBLIC
    OpenMP::OpenMP_CXX
)

if(MSVC)
    target_compile_options(sonar_core PRIVATE /W4 /O2)
else()
    target_compile_options(sonar_core PRIVATE -Wall -Wextra -O3 -march=native)
endif()

if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.12.0
    )
    FetchContent_MakeAvailable(pybind11)

    pybind11_add_module(sonar_py
        python/bindings.cpp
    )

    target_link_libraries(sonar_py PRIVATE sonar_core)

    set_target_properties(sonar_py PROPERTIES
        OUTPUT_NAME sonar
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
    )
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
