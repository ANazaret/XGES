//
// Created by Achille Nazaret on 11/4/23.
//

#include "Operators.h"


Insert::Insert(int x, int y, FlatSet T, double score, FlatSet effective_parents)
    : x(x), y(y), T(std::move(T)), score(score),
      effective_parents(std::move(effective_parents)) {}

std::ostream &operator<<(std::ostream &os, const Insert &obj) {
    os << "Insert: " << obj.x << " → " << obj.y << ", T = {";
    for (auto t: obj.T) { os << t << ", "; }
    os << "}, score = " << obj.score << ", effective_parents = {";
    for (auto p: obj.effective_parents) { os << p << ", "; }
    os << "}";
    return os;
}


Delete::Delete(int x, int y, FlatSet C, double score, FlatSet effective_parents,
               bool directed)
    : x(x), y(y), C(std::move(C)), score(score),
      effective_parents(std::move(effective_parents)), directed(directed) {}

std::ostream &operator<<(std::ostream &os, const Delete &obj) {
    os << "Delete: " << obj.x << (obj.directed ? " → " : " - ") << obj.y << ", O = {";
    for (auto t: obj.C) { os << t << ", "; }
    os << "}, score = " << obj.score;
    os << ", effective_parents = {";
    for (auto p: obj.effective_parents) { os << p << ", "; }
    os << "}";
    return os;
}

Reverse::Reverse(Insert insert, double score, FlatSet parents_x)
    : insert(std::move(insert)), score(score), parents_x(std::move(parents_x)) {}

Reverse::Reverse(int x, int y, const FlatSet &T, double score,
                 const FlatSet &effective_parents, FlatSet parents_x)
    : insert(x, y, T, 0, effective_parents), score(score),
      parents_x(std::move(parents_x)) {}

std::ostream &operator<<(std::ostream &os, const Reverse &obj) {
    os << "Reverse: " << obj.insert.x << " ← " << obj.insert.y << ", T = {";
    for (auto t: obj.insert.T) { os << t << ", "; }
    os << "}, score = " << obj.score << ", effective_parents = {";
    for (auto p: obj.insert.effective_parents) { os << p << ", "; }
    os << "}";
    return os;
}
