cmake_minimum_required(VERSION 3.9)
project(molecular_dynamics_poc)

set(CMAKE_CXX_STANDARD 11)

add_executable(molecular_dynamics_poc src/main.cc
        src/vectors.cc src/vectors.h
        src/parameters.cc src/parameters.h src/atoms.cc src/atoms.h)
