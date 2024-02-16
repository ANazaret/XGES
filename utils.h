//
// Created by Achille Nazaret on 2/16/24.
//

#pragma once

#include <chrono>
using namespace std::chrono;

inline double measure_time(const high_resolution_clock::time_point start_time) {
    return duration_cast<duration<double>>(high_resolution_clock::now() - start_time)
            .count();
}