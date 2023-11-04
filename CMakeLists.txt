cmake_minimum_required(VERSION 3.26)
project(ges_cpp)

set(CMAKE_CXX_STANDARD 23)

add_executable(ges_cpp main.cpp
        PDAG.cpp
        PDAG.h)
