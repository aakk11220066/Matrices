cmake_minimum_required(VERSION 3.13)
project(Matrices)

set(CMAKE_CXX_STANDARD 11)

add_executable(Matrices
        Auxilaries.h
        Complex.h
        MtmExceptions.h
        MtmMat.h
        MtmMatSq.h
        MtmMatTriag.h
        MtmVec.h
        main.cpp)