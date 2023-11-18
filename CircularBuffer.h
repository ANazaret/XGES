//
// Created by Achille Nazaret on 11/17/23.
//

#pragma once

#include <memory>

template<class T>
class CircularBuffer {
private:
    std::unique_ptr<T[]> buffer;
    int front_idx = 0;
    int back_idx = 0;
    int max_size;

public:
    CircularBuffer<T>(int max_size) : buffer(std::unique_ptr<T[]>(new T[max_size])), max_size(max_size){};

    void push_back(T item) {
        buffer[back_idx] = item;
        back_idx = (back_idx + 1) % max_size;
    }

    T pop_back() {
        if (empty()) throw std::runtime_error("CircularBuffer is empty");

        // move back_idx back
        back_idx = (back_idx - 1) % max_size;
        return buffer[back_idx];
    }

    T pop_front() {
        if (empty()) throw std::runtime_error("CircularBuffer is empty");
        T item = buffer[front_idx];
        front_idx = (front_idx + 1) % max_size;
        return item;
    }

    bool empty() { return front_idx == back_idx; }

    int size() {
        if (front_idx >= back_idx) return front_idx - back_idx;
        else
            return max_size - back_idx + front_idx;
    }

    void clear() {
        front_idx = 0;
        back_idx = 0;
    }
};