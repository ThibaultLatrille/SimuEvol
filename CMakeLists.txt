cmake_minimum_required(VERSION 3.7)
project(SimuEvol)

set(CMAKE_CXX_STANDARD 14)

# Compilation options
option(COVERAGE_MODE "For coverage mode using g++ " OFF) #OFF by default
if(COVERAGE_MODE)
    set(CMAKE_CXX_FLAGS "-O0 --std=c++11 -fprofile-arcs -ftest-coverage") # coverage mode
else()
    set(CMAKE_CXX_FLAGS "-O0 --std=c++11 -Wall -Wextra -Wno-unused-parameter -Wfatal-errors") # normal mode
endif(COVERAGE_MODE)
unset(COVERAGE_MODE CACHE)

include_directories("src")
include_directories("src/lib")
include_directories("utils")

add_executable(SimuEvol src/simu_evol.cpp )
add_executable(SimuPoly src/simu_poly.cpp)
add_executable(SimuRelax src/simu_relax.cpp)