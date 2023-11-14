//
// Created by Achille Nazaret on 11/6/23.
//

#pragma once

#include <queue>
#include <set>

enum class EdgeType { DIRECTED = 1, UNDIRECTED = 2 };

struct Edge {
    int x, y;
    EdgeType type;

    bool operator<(const Edge &other) const {
        // TODO: change to make undirected edges equal to their reverse
        if (x < other.x) return true;
        if (x > other.x) return false;
        if (y < other.y) return true;
        if (y > other.y) return false;
        return type < other.type;
    }
};

class EdgeQueueSet {
private:
    std::queue<Edge> edges_queue;
    std::set<Edge> edges_set;

public:
    void push_directed(int x, int y);

    void push_undirected(int x, int y);

    Edge pop();

    bool empty() const;
};
