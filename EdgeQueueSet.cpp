//
// Created by Achille Nazaret on 11/6/23.
//

#include "EdgeQueueSet.h"
#include <iostream>

Edge::Edge(int x, int y, EdgeType type) {
    if (x > y) {
        std::swap(x, y);
        if (type == EdgeType::DIRECTED_TO_X) type = EdgeType::DIRECTED_TO_Y;
        else if (type == EdgeType::DIRECTED_TO_Y)
            type = EdgeType::DIRECTED_TO_X;
    }
    this->x = x;
    this->y = y;
    this->type = type;
}

void EdgeQueueSet::push_directed(int x, int y) {
    Edge e{x, y, EdgeType::DIRECTED_TO_Y};
    if (!edges_set.contains(e)) {
        edges_queue.push(e);
        edges_set.insert(e);
    }
}

void EdgeQueueSet::push_undirected(int x, int y) {
    if (x > y) std::swap(x, y);
    Edge e{x, y, EdgeType::UNDIRECTED};
    if (!edges_set.contains(e)) {
        edges_queue.push(e);
        edges_set.insert(e);
    }
}

Edge EdgeQueueSet::pop() {
    Edge e = edges_queue.front();
    edges_queue.pop();
    edges_set.erase(e);
    return e;
}

bool EdgeQueueSet::empty() const { return edges_queue.empty(); }


EdgeModification::EdgeModification(int x, int y, EdgeType old_type, EdgeType new_type)
    : x(x), y(y), old_type(old_type), new_type(new_type) {}

bool EdgeModification::is_reverse() const {
    return (old_type == EdgeType::DIRECTED_TO_X && new_type == EdgeType::DIRECTED_TO_Y) ||
           (old_type == EdgeType::DIRECTED_TO_Y && new_type == EdgeType::DIRECTED_TO_X);
}

bool EdgeModification::is_new_directed() const {
    return new_type == EdgeType::DIRECTED_TO_X || new_type == EdgeType::DIRECTED_TO_Y;
}

bool EdgeModification::is_old_directed() const {
    return old_type == EdgeType::DIRECTED_TO_X || old_type == EdgeType::DIRECTED_TO_Y;
}

bool EdgeModification::is_new_undirected() const {
    return new_type == EdgeType::UNDIRECTED;
}

bool EdgeModification::is_old_undirected() const {
    return old_type == EdgeType::UNDIRECTED;
}

int EdgeModification::get_new_target() const {
    if (new_type == EdgeType::DIRECTED_TO_X) return x;
    if (new_type == EdgeType::DIRECTED_TO_Y) return y;
    throw std::runtime_error("New edge is not directed");
}

int EdgeModification::get_new_source() const {
    if (new_type == EdgeType::DIRECTED_TO_X) return y;
    if (new_type == EdgeType::DIRECTED_TO_Y) return x;
    throw std::runtime_error("New edge is not directed");
}

int EdgeModification::get_old_target() const {
    if (old_type == EdgeType::DIRECTED_TO_X) return x;
    if (old_type == EdgeType::DIRECTED_TO_Y) return y;
    throw std::runtime_error("Old edge is not directed");
}

int EdgeModification::get_old_source() const {
    if (old_type == EdgeType::DIRECTED_TO_X) return y;
    if (old_type == EdgeType::DIRECTED_TO_Y) return x;
    throw std::runtime_error("Old edge is not directed");
}
int EdgeModification::get_modification_id() const {
    if (old_type == EdgeType::NONE && is_new_undirected()) { return 1; }
    if (old_type == EdgeType::NONE && is_new_directed()) { return 2; }
    if (is_old_undirected() && new_type == EdgeType::NONE) { return 3; }
    if (is_old_undirected() && is_new_directed()) { return 4; }
    if (is_old_directed() && new_type == EdgeType::NONE) { return 5; }
    if (is_old_directed() && is_new_undirected()) { return 6; }
    if (is_reverse()) { return 7; }
    throw std::runtime_error("Invalid modification");
}


void EdgeModificationsMap::update_edge_directed(int x, int y, EdgeType old_type) {
    EdgeType new_type;
    if (x > y) {
        std::swap(x, y);
        new_type = EdgeType::DIRECTED_TO_X;
        if (old_type == EdgeType::DIRECTED_TO_X) old_type = EdgeType::DIRECTED_TO_Y;
    } else {
        new_type = EdgeType::DIRECTED_TO_Y;
    }
    update_edge_modification(x, y, old_type, new_type);
}

void EdgeModificationsMap::update_edge_undirected(int x, int y, EdgeType old_type) {
    if (x > y) {
        std::swap(x, y);
        if (old_type == EdgeType::DIRECTED_TO_X) old_type = EdgeType::DIRECTED_TO_Y;
        else if (old_type == EdgeType::DIRECTED_TO_Y)
            old_type = EdgeType::DIRECTED_TO_X;
    }
    update_edge_modification(x, y, old_type, EdgeType::UNDIRECTED);
}

void EdgeModificationsMap::update_edge_none(int x, int y, EdgeType old_type) {
    if (x > y) {
        std::swap(x, y);
        if (old_type == EdgeType::DIRECTED_TO_X) old_type = EdgeType::DIRECTED_TO_Y;
        else if (old_type == EdgeType::DIRECTED_TO_Y)
            old_type = EdgeType::DIRECTED_TO_X;
    }
    update_edge_modification(x, y, old_type, EdgeType::NONE);
}


void EdgeModificationsMap::update_edge_modification(int small, int big, EdgeType old_type,
                                                    EdgeType new_type) {
    std::pair<int, int> edge_key(small, big);

    if (const auto edge_modification = edge_modifications.find(edge_key);
        edge_modification != edge_modifications.end()) {
        // the edge was already modified

        auto oldest_type = edge_modification->second.old_type;
        if (oldest_type == new_type) {
            // we delete the edge as it is not modified anymore
            edge_modifications.erase(edge_key);
        } else {
            // we update the type of the edge modification
            edge_modification->second.new_type = new_type;
        }
    } else {
        // we add the edge to the list of modified edges
        edge_modifications.emplace(edge_key,
                                   EdgeModification(small, big, old_type, new_type));
    }
}

void EdgeModificationsMap::clear() { edge_modifications.clear(); }

std::map<std::pair<int, int>, EdgeModification>::iterator EdgeModificationsMap::begin() {
    return edge_modifications.begin();
}

std::map<std::pair<int, int>, EdgeModification>::iterator EdgeModificationsMap::end() {
    return edge_modifications.end();
}

std::ostream &operator<<(std::ostream &os, const EdgeType &edge_type) {
    switch (edge_type) {
        case EdgeType::DIRECTED_TO_X:
            os << "←";
            break;
        case EdgeType::DIRECTED_TO_Y:
            os << "→";
            break;
        case EdgeType::UNDIRECTED:
            os << "-";
            break;
        case EdgeType::NONE:
            os << " ";
            break;
    }
    return os;
}

// print function for EdgeModification
std::ostream &operator<<(std::ostream &os, const EdgeModification &edge_modification) {
    os << edge_modification.x << " " << edge_modification.old_type << " "
       << edge_modification.y << " becomes " << edge_modification.x << " "
       << edge_modification.new_type << " " << edge_modification.y;
    return os;
}