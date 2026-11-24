# Python extension: nanobind bindings over the C++ core at the repository root
# (built from the same sources, with the same strict floating-point flags
# (see ../cpp/cmake/LexFloatingPoint.cmake), so results are bit-identical to the C++ library).
cmake_minimum_required(VERSION 3.21)
project(lex1dx_python LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(LEX1DX_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH "lex1dx repository root")
if(NOT EXISTS "${LEX1DX_ROOT}/CMakeLists.txt" OR NOT EXISTS "${LEX1DX_ROOT}/cpp/include/lex1dx/lex1dx.hpp")
  message(FATAL_ERROR "C++ core not found at ${LEX1DX_ROOT}. Build from a full repository checkout.")
endif()
set(LEX1DX_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LEX1DX_BUILD_CLI OFF CACHE BOOL "" FORCE)
add_subdirectory("${LEX1DX_ROOT}" lex1dx_core EXCLUDE_FROM_ALL)

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(_core NB_STATIC src/_core.cpp)
target_link_libraries(_core PRIVATE lex1dx::lex1dx)
include("${LEX1DX_ROOT}/cpp/cmake/LexFloatingPoint.cmake")
lex_strict_fp(_core)
install(TARGETS _core LIBRARY DESTINATION lex1dx)
