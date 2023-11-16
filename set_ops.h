//
// Created by Achille Nazaret on 11/4/23.
//

#pragma once

#include <iostream>
#include <set>

bool is_subset(const std::set<int> &a, const std::set<int> &b);

bool have_overlap(const std::set<int> &a, const std::set<int> &b);

void intersect_in_place(std::set<int> &a, const std::set<int> &b);

bool equal_union(const std::set<int> &a, const std::set<int> &b1, const std::set<int> &b2);

// printer of std::set<int>
std::ostream &operator<<(std::ostream &os, const std::set<int> &obj);
