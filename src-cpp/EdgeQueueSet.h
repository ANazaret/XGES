//
// Created by Achille Nazaret on 11/6/23.
//

#pragma once
#include <map>
#include <queue>
#include <set>

enum class EdgeType { NONE = 0, UNDIRECTED = 1, DIRECTED_TO_Y = 2, DIRECTED_TO_X = 3 };

class Edge {
    int x, y;
    EdgeType type;

public:
    bool operator<(const Edge &other) const {
        if (x < other.x) return true;
        if (x > other.x) return false;
        if (y < other.y) return true;
        if (y > other.y) return false;
        return type < other.type;
    }

    Edge(int x, int y, EdgeType type);

    bool is_directed() const {
        return (type == EdgeType::DIRECTED_TO_X || type == EdgeType::DIRECTED_TO_Y);
    }
    int get_source() const {
        if (type == EdgeType::DIRECTED_TO_X) return y;
        if (type == EdgeType::DIRECTED_TO_Y) return x;
        throw std::runtime_error("Edge is not directed");
    }
    int get_target() const {
        if (type == EdgeType::DIRECTED_TO_X) return x;
        if (type == EdgeType::DIRECTED_TO_Y) return y;
        throw std::runtime_error("Edge is not directed");
    }
    int get_x() const { return x; }
    int get_y() const { return y; }
};

class EdgeQueueSet {
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
    int x, y;
    EdgeType old_type, new_type;

    EdgeModification(int x, int y, EdgeType old_type, EdgeType new_type);

    bool is_reverse() const;

    bool is_new_directed() const;
    bool is_old_directed() const;

    bool is_new_undirected() const;
    bool is_old_undirected() const;

    int get_new_target() const;

    int get_new_source() const;
    int get_old_target() const;
    int get_old_source() const;

    int get_modification_id() const;
};


class EdgeModificationsMap {
public:
    void update_edge_directed(int x, int y, EdgeType old_type);
    void update_edge_undirected(int x, int y, EdgeType old_type);
    void update_edge_none(int x, int y, EdgeType old_type);
    void clear();

    std::map<std::pair<int, int>, EdgeModification>::iterator begin();

    std::map<std::pair<int, int>, EdgeModification>::iterator end();


private:
    void update_edge_modification(int small, int big, EdgeType old_type,
                                  EdgeType new_type);
    std::map<std::pair<int, int>, EdgeModification> edge_modifications;
};

std::ostream &operator<<(std::ostream &os, const EdgeModification &edge_modification);
