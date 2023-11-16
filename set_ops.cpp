//
// Created by Achille Nazaret on 11/4/23.
//

#include "set_ops.h"

bool is_subset(const std::set<int> &a, const std::set<int> &b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

bool have_overlap(const std::set<int> &a, const std::set<int> &b) {
    if (a.empty() || b.empty()) { return false; }
    if (a.size() > b.size()) { return have_overlap(b, a); }
    for (auto &i: a) {
        if (b.find(i) != b.end()) { return true; }
    }
    return false;
}

/**
 * Intersect a and b in place.
 *
 * From https://stackoverflow.com/questions/1773526/in-place-c-set-intersection
 * @param a
 * @param b
 */
void intersect_in_place(std::set<int> &a, const std::set<int> &b) {
    std::set<int>::iterator it1 = a.begin();
    std::set<int>::iterator it2 = b.begin();
    while ((it1 != a.end()) && (it2 != b.end())) {
        if (*it1 < *it2) {
            // a.erase(it1++);
            it1 = a.erase(it1);
        } else if (*it2 < *it1) {
            ++it2;
        } else {// *it1 == *it2
            ++it1;
            ++it2;
        }
    }
    a.erase(it1, a.end());
}

/**
 * Check if a equals the union of b1 and b2.
 * Does not perform an explicit union for efficiency [i hope it's more efficient].
 *
 * @param a
 * @param b1
 * @param b2
 * @return
 */
bool equal_union(const std::set<int> &a, const std::set<int> &b1, const std::set<int> &b2) {
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

std::ostream &operator<<(std::ostream &os, const std::set<int> &obj) {
    os << "{";
    for (auto &i: obj) { os << i << ", "; }
    os << "}";
    return os;
}