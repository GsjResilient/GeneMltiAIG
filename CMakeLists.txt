cmake_minimum_required(VERSION 3.16)

project(GeneMltiAIG VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ec_expand
    src/aag.cpp
    src/filters.cpp
    src/main.cpp
    src/validation.cpp
)

target_include_directories(ec_expand PRIVATE include)

if (MSVC)
    target_compile_options(ec_expand PRIVATE /W4 /permissive-)
else()
    target_compile_options(ec_expand PRIVATE -Wall -Wextra -Wpedantic)
endif()
