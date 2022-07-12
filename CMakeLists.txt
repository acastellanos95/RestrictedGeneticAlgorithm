cmake_minimum_required(VERSION 3.22)
project(RestrictedGeneticAlgorithm)

set(CMAKE_CXX_STANDARD 17)

add_executable(RestrictedGeneticAlgorithm main.cpp Utils.h Individual.h lib/RandomGA.h)
