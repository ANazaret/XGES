cmake_minimum_required(VERSION 3.26)
project(ges_cpp)

set(CMAKE_CXX_STANDARD 23)
set(EIGEN3_INCLUDE_DIR "eigen-3.4.0")
include_directories(${EIGEN3_INCLUDE_DIR})
add_subdirectory(cnpy)

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
)

# Include library headers
# Include the cnpy directory
target_include_directories(ges_cpp PRIVATE cnpy)
target_link_libraries(ges_cpp PRIVATE cnpy)


find_package(Boost REQUIRED)
target_link_libraries(ges_cpp PRIVATE Boost::boost)

find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(ges_cpp PRIVATE ${BLAS_LIBRARIES})
