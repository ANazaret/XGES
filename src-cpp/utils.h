//
// Created by Achille Nazaret on 2/16/24.
//

#pragma once

#include <chrono>

inline double
measure_time(const std::chrono::high_resolution_clock::time_point start_time) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
                   std::chrono::high_resolution_clock::now() - start_time)
            .count();
}

// define hash function for pairs of ints
template<>
struct std::hash<std::pair<int, int>> {
    size_t operator()(const std::pair<int, int> &p) const noexcept {
        size_t h1 = std::hash<int>()(p.first);
        h1 ^= std::hash<int>()(p.second) + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
        return h1;
    }
};// namespace std
