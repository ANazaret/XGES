//
// Created by Achille Nazaret on 11/4/23.
//

#include "Insert.h"

Insert::Insert(int x, int y, const std::set<int> &T, double score, std::set<int> effective_parents)
    : x(x), y(y), T(T), score(score), effective_parents(effective_parents) {}

std::ostream &operator<<(std::ostream &os, const Insert &obj) {
    os << "Insert: x = " << obj.x << ", y = " << obj.y << ", T = {";
    for (auto t: obj.T) { os << t << ", "; }
    os << "}, score = " << obj.score << ", effective_parents = {";
    for (auto p: obj.effective_parents) { os << p << ", "; }
    os << "}";
    return os;
}