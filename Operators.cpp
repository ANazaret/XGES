//
// Created by Achille Nazaret on 11/4/23.
//

#include "Operators.h"

Insert::Insert(int x, int y, const std::set<int> &T, double score, std::set<int> effective_parents)
    : x(x), y(y), T(T), score(score), effective_parents(effective_parents) {}

std::ostream &operator<<(std::ostream &os, const Insert &obj) {
    os << "Insert: " << obj.x << " → " << obj.y << ", T = {";
    for (auto t: obj.T) { os << t << ", "; }
    os << "}, score = " << obj.score << ", effective_parents = {";
    for (auto p: obj.effective_parents) { os << p << ", "; }
    os << "}";
    return os;
}


Delete::Delete(int x, int y, const std::set<int> &O, double score, std::set<int> effective_parents)
    : x(x), y(y), O(O), score(score) {}

std::ostream &operator<<(std::ostream &os, const Delete &obj) {
    os << "Delete: " << obj.x << "-?-" << obj.y << ", O = {";
    for (auto t: obj.O) { os << t << ", "; }
    os << "}, score = " << obj.score;
    return os;
}

Reverse::Reverse(Insert insert, double score) : insert(insert), score(score) {}

Reverse::Reverse(int x, int y, const std::set<int> &T, double score, std::set<int> effective_parents)
    : insert(x, y, T, 0, effective_parents), score(score) {}

std::ostream &operator<<(std::ostream &os, const Reverse &obj) {
    os << "Reverse: " << obj.insert.x << " ← " << obj.insert.y << ", T = {";
    for (auto t: obj.insert.T) { os << t << ", "; }
    os << "}, score = " << obj.score << ", effective_parents = {";
    for (auto p: obj.insert.effective_parents) { os << p << ", "; }
    os << "}";
    return os;
}