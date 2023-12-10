//
// Created by Achille Nazaret on 11/17/23.
//

#pragma once

#include <memory>

template<class T>
class CircularBuffer {
private:
    std::unique_ptr<T[]> buffer;
    int front_idx = 0;// index of the first element
    int back_idx = 0; // index of the last element
    int max_size;

public:
    CircularBuffer<T>(int max_size) : buffer(std::unique_ptr<T[]>(new T[max_size])), max_size(max_size){};

    // copy assignment
    CircularBuffer<T> &operator=(const CircularBuffer<T> &other) {
        buffer = std::unique_ptr<T[]>(new T[other.max_size]);
        max_size = other.max_size;
        front_idx = other.front_idx;
        back_idx = other.back_idx;
        if (front_idx < back_idx)
            std::copy(other.buffer.get() + front_idx, other.buffer.get() + back_idx, buffer.get());
        else if (front_idx > back_idx) {
            std::copy(other.buffer.get() + front_idx, other.buffer.get() + max_size, buffer.get());
            std::copy(other.buffer.get(), other.buffer.get() + back_idx, buffer.get() + max_size - front_idx);
        }
        return *this;
    }

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
        if (front_idx <= back_idx) return back_idx - front_idx;
        else
            return max_size + back_idx - front_idx;
    }

    void clear() {
        front_idx = 0;
        back_idx = 0;
    }

    CircularBuffer<T>(const CircularBuffer<T> &other)
        : buffer(std::unique_ptr<T[]>(new T[other.max_size])), max_size(other.max_size) {
        front_idx = other.front_idx;
        back_idx = other.back_idx;
        std::copy(other.buffer.get(), other.buffer.get() + max_size, buffer.get());
    }
};