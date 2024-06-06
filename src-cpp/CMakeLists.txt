set(BINARY ${CMAKE_PROJECT_NAME})

add_executable(
        ${BINARY} main.cpp
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
        utils.h
)

target_link_libraries(${BINARY} PRIVATE cnpy)
target_link_libraries(${BINARY} PRIVATE spdlog::spdlog)

find_package(Boost REQUIRED)
target_link_libraries(${BINARY} PRIVATE Boost::boost)

find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(${BINARY} PRIVATE ${BLAS_LIBRARIES})
