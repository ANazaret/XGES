//
// Created by Achille Nazaret on 11/6/23.
//

#pragma once

#include <queue>
#include <set>

enum class EdgeType { UNDIRECTED = 2, NONE = 3, DIRECTED_TO_X = 4, DIRECTED_TO_Y = 5 };

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


class EdgeModification {
public:
    int x, y;// always x < y
    EdgeType old_type, new_type;

    EdgeModification(int x, int y, EdgeType old_type, EdgeType new_type);

    bool is_now_reverse() const;

    bool is_now_directed() const;

    bool is_now_undirected() const;

    int get_target() const;

    int get_source() const;
};

#include <map>

class EdgeModificationsMap {
public:
    void update_edge_directed(int x, int y, EdgeType old_type);
    void update_edge_undirected(int x, int y, EdgeType old_type);
    void update_edge_none(int x, int y, EdgeType old_type);

    std::vector<EdgeModification> get_edge_modifications() const;

    void clear();

    std::map<std::pair<int, int>, EdgeModification>::iterator begin();
    
    std::map<std::pair<int, int>, EdgeModification>::iterator end();


private:
    void update_edge_modification(int small, int big, EdgeType old_type, EdgeType new_type);
    std::map<std::pair<int, int>, EdgeModification> edge_modifications;
};