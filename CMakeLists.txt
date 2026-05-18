cmake_minimum_required(VERSION 3.20)

project(flash_cubical LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core
    cpp/bindings.cpp
    cpp/compute_2d_all.cpp
    cpp/zero_lookup/runtime/mask_build.cpp
    cpp/zero_lookup/runtime/zero_table.cpp
)

target_include_directories(_core PRIVATE cpp)

if (MSVC)
    target_compile_options(_core PRIVATE /O2 /W4 /permissive-)
else()
    target_compile_options(_core PRIVATE -O3 -Wall -Wextra)
endif()

install(TARGETS _core DESTINATION flash_cubical)
