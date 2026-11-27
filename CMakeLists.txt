cmake_minimum_required(VERSION 3.21)
project(lex1dx VERSION 0.1.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cpp/cmake/LexFloatingPoint.cmake)

option(LEX1DX_BUILD_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
option(LEX1DX_BUILD_CLI "Build the lex1dx_dump tool" ${PROJECT_IS_TOP_LEVEL})

add_library(lex1dx
  cpp/src/atoms.cpp
  cpp/src/collapse.cpp
  cpp/src/refine.cpp
  cpp/src/exact.cpp
  cpp/src/sort.cpp
  cpp/src/lex.cpp
  cpp/src/lex1dx.cpp
  third_party/kmeans1d/interval_sum.cpp
  third_party/kmeans1d/kmeans_wilber.cpp
)
add_library(lex1dx::lex1dx ALIAS lex1dx)
set_target_properties(lex1dx PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(lex1dx
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${PROJECT_SOURCE_DIR}/cpp/src
    ${PROJECT_SOURCE_DIR}/third_party/kmeans1d
)
lex_strict_fp(lex1dx)
if(MSVC)
  target_compile_options(lex1dx PRIVATE /W4)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(lex1dx PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(LEX1DX_BUILD_CLI)
  add_executable(lex1dx_dump cpp/tools/lex1dx_dump.cpp)
  target_link_libraries(lex1dx_dump PRIVATE lex1dx)
  lex_strict_fp(lex1dx_dump)
endif()

if(LEX1DX_BUILD_TESTS)
  enable_testing()
  foreach(t edge quality agree)
    add_executable(test_${t} cpp/tests/test_${t}.cpp)
    target_link_libraries(test_${t} PRIVATE lex1dx)
    target_include_directories(test_${t} PRIVATE ${PROJECT_SOURCE_DIR}/cpp/src ${PROJECT_SOURCE_DIR}/cpp/include)
    lex_strict_fp(test_${t})
    add_test(NAME ${t} COMMAND test_${t})
  endforeach()
  # exact.cpp is part of the library; quality/edge tests may also call lex::ckmeans_exact
endif()
