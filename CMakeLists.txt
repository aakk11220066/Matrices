cmake_minimum_required(VERSION 3.13)
project(Matrices)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -Wall -Werror -pedantic-errors")

add_executable(Matrices
        Auxilaries.h
        Complex.h
        MtmExceptions.h
        MtmMat.h
        MtmMatSq.h
        MtmMatTriag.h
        MtmVec.h
        RootVector.h
        main.cpp)