cmake_minimum_required(VERSION 3.16)

project(GeneMltiAIG VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ec_expand
    src/aag.cpp
    src/filters.cpp
    src/main.cpp
    src/sub_aig.cpp
    src/validation.cpp
)

add_executable(aig_sub_enum
    src/aag.cpp
    src/aig_sub_enum.cpp
    src/sub_aig.cpp
)

target_include_directories(ec_expand PRIVATE include)
target_include_directories(aig_sub_enum PRIVATE include)

if (MSVC)
    target_compile_options(ec_expand PRIVATE /W4 /permissive-)
    target_compile_options(aig_sub_enum PRIVATE /W4 /permissive-)
else()
    target_compile_options(ec_expand PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(aig_sub_enum PRIVATE -Wall -Wextra -Wpedantic)
endif()
