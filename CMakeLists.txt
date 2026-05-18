cmake_minimum_required(VERSION 3.20)

project(flash_cubical LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FLASH_CUBICAL_LOOKUP_DATA_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/src/flash_cubical/data"
    CACHE PATH "Directory containing packaged zero lookup table files")

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core
    cpp/bindings.cpp
    cpp/compute_2d_all.cpp
    cpp/zero_lookup/runtime/mask_build.cpp
    cpp/zero_lookup/runtime/zero_table.cpp
    cpp/zero_lookup/runtime/zero_table_file.cpp
)

target_include_directories(_core PRIVATE cpp src .)

if (MSVC)
    target_compile_options(_core PRIVATE /O2 /W4 /permissive-)
else()
    target_compile_options(_core PRIVATE -O3 -Wall -Wextra)
endif()

install(TARGETS _core DESTINATION flash_cubical)
install(FILES
    "${FLASH_CUBICAL_LOOKUP_DATA_DIR}/zero_lookup_2d.bin"
    "${FLASH_CUBICAL_LOOKUP_DATA_DIR}/zero_lookup_3d.bin"
    DESTINATION flash_cubical/data)
