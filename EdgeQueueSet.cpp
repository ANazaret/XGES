//
// Created by Achille Nazaret on 11/6/23.
//

#include "EdgeQueueSet.h"

void EdgeQueueSet::push_directed(int x, int y) {
    Edge e{x, y, EdgeType::DIRECTED_TO_Y};
    if (edges_set.find(e) == edges_set.end()) {
        edges_queue.push(e);
        edges_set.insert(e);
    }
}

void EdgeQueueSet::push_undirected(int x, int y) {
    if (x > y) std::swap(x, y);
    Edge e{x, y, EdgeType::UNDIRECTED};
    if (edges_set.find(e) == edges_set.end()) {
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

bool EdgeModification::is_now_reverse() const {
    return (old_type == EdgeType::DIRECTED_TO_X && new_type == EdgeType::DIRECTED_TO_Y) ||
           (old_type == EdgeType::DIRECTED_TO_Y && new_type == EdgeType::DIRECTED_TO_X);
}

bool EdgeModification::is_now_directed() const {
    return new_type == EdgeType::DIRECTED_TO_X || new_type == EdgeType::DIRECTED_TO_Y;
}

bool EdgeModification::is_now_undirected() const { return new_type == EdgeType::UNDIRECTED; }

int EdgeModification::get_target() const {
    if (new_type == EdgeType::DIRECTED_TO_X) return x;
    if (new_type == EdgeType::DIRECTED_TO_Y) return y;
    return -1;
}

int EdgeModification::get_source() const {
    if (new_type == EdgeType::DIRECTED_TO_X) return y;
    if (new_type == EdgeType::DIRECTED_TO_Y) return x;
    return -1;
}


void EdgeModificationsMap::update_edge_directed(int x, int y, EdgeType old_type) {
    EdgeType new_type;
    if (x > y) {
        std::swap(x, y);
        new_type = EdgeType::DIRECTED_TO_X;
        assert(old_type != EdgeType::DIRECTED_TO_Y);
        if (old_type == EdgeType::DIRECTED_TO_X) old_type = EdgeType::DIRECTED_TO_Y;
    } else {
        assert(old_type != EdgeType::DIRECTED_TO_Y);
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
            assert(edge_modification->second.new_type == old_type);
            edge_modification->second.new_type = new_type;
        }
    } else {
        // we add the edge to the list of modified edges
        edge_modifications.emplace(edge_key, EdgeModification(small, big, old_type, new_type));
    }
}

void EdgeModificationsMap::clear() { edge_modifications.clear(); }

std::map<std::pair<int, int>, EdgeModification>::iterator EdgeModificationsMap::begin() {
    return edge_modifications.begin();
}

std::map<std::pair<int, int>, EdgeModification>::iterator EdgeModificationsMap::end() {
    return edge_modifications.end();
}