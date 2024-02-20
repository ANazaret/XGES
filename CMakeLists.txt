cmake_minimum_required(VERSION 3.26)
project(ges_cpp)

set(CMAKE_CXX_STANDARD 23)
include_directories(Eigen)

add_executable(
        ges_cpp main.cpp
        PDAG.cpp
        PDAG.h
        ScorerInterface.h
        Operators.h
        Operators.cpp
        set_ops.h
        set_ops.cpp
        EdgeQueueSet.cpp
        EdgeQueueSet.h
        XGES.cpp
        XGES.h
        BICScorer.cpp
        BICScorer.h
        CircularBuffer.h
        dependencies/cxxopts.hpp
        utils.h
)

add_subdirectory(cnpy)
target_include_directories(ges_cpp PRIVATE cnpy)
target_link_libraries(ges_cpp PRIVATE cnpy)

add_subdirectory(spdlog)
target_link_libraries(ges_cpp PRIVATE spdlog::spdlog)

find_package(Boost REQUIRED)
target_link_libraries(ges_cpp PRIVATE Boost::boost)

find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(ges_cpp PRIVATE ${BLAS_LIBRARIES})
