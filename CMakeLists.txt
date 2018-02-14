cmake_minimum_required(VERSION 3.7)
project(SimuEvol)

set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES Eigen/Dense docopt.cpp/docopt.h codon.hpp main.cpp )

add_executable(SimuEvol ${SOURCE_FILES})