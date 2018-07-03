cmake_minimum_required(VERSION 3.7)
project(SimuEvol)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES_EVOL src/Eigen/Dense src/docopt.cpp/docopt.h src/codon.hpp src/simu_evol.cpp )
add_executable(SimuEvol ${SOURCE_FILES_EVOL})

set(SOURCE_FILES_POLY src/docopt.cpp/docopt.h src/codon.hpp src/simu_poly.cpp )
add_executable(SimuPoly ${SOURCE_FILES_POLY})

set(SOURCE_FILES_RELAX src/docopt.cpp/docopt.h src/simu_relax.cpp )
add_executable(SimuRelax ${SOURCE_FILES_RELAX})