//
// Created by Achille Nazaret on 11/6/23.
//

#include "EdgeQueueSet.h"

void EdgeQueueSet::push_directed(int x, int y) {
    Edge e{x, y, EdgeType::DIRECTED};
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

bool EdgeQueueSet::empty() const {
    return edges_queue.empty();
}
