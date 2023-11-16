//
// Created by Achille Nazaret on 11/4/23.
//

#include "Operators.h"

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


Delete::Delete(int x, int y, const std::set<int> &O, double score, std::set<int> effective_parents)
    : x(x), y(y), O(O), score(score) {}

std::ostream &operator<<(std::ostream &os, const Delete &obj) {
    os << "Delete: x = " << obj.x << ", y = " << obj.y << ", H = {";
    for (auto t: obj.H) { os << t << ", "; }
    os << "}, score = " << obj.score;
    return os;
}