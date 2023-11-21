//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <iostream>
#include <set>

typedef boost::container::flat_set<int, std::less<int>, boost::container::small_vector<int, 10>> FlatSet;


template<typename T>
bool is_subset(const T &a, const T &b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

template<typename T>
bool have_overlap(const T &a, const T &b) {
    if (a.empty() || b.empty()) { return false; }
    if (a.size() > b.size()) { return have_overlap(b, a); }
    for (auto &i: a) {
        if (b.find(i) != b.end()) { return true; }
    }
    return false;
}


template<typename T>
bool equal_union(const T &a, const T &b1, const T &b2) {
    // Iterate over a and check that each element is in b1 or b2
    auto itb1 = b1.begin();
    auto itb2 = b2.begin();
    auto ita = a.begin();

    while (ita != a.end()) {
        if (itb1 != b1.end() && *itb1 == *ita) {
            ++itb1;
            ++ita;
        } else if (itb2 != b2.end() && *itb2 == *ita) {
            ++itb2;
            ++ita;
        } else {
            return false;
        }
    }
    if (itb1 != b1.end() || itb2 != b2.end()) { return false; }
    return true;
}

void union_with_single_element(const FlatSet &a, int b, FlatSet &out);