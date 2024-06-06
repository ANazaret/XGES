from enum import Enum
from collections import deque


class EdgeType(Enum):
    NONE = 0
    UNDIRECTED = 1
    DIRECTED_TO_Y = 2
    DIRECTED_TO_X = 3

    @staticmethod
    def to_string(edge_type):
        if edge_type == EdgeType.NONE:
            return "*"
        if edge_type == EdgeType.UNDIRECTED:
            return "-"
        if edge_type == EdgeType.DIRECTED_TO_Y:
            return "→"
        if edge_type == EdgeType.DIRECTED_TO_X:
            return "←"


class Edge:
    def __init__(self, x, y, type):
        self.x = x
        self.y = y
        self.type = type

    def __lt__(self, other):
        if self.x < other.x:
            return True
        if self.x > other.x:
            return False
        if self.y < other.y:
            return True
        if self.y > other.y:
            return False
        return self.type < other.type

    def is_directed(self):
        return self.type == EdgeType.DIRECTED_TO_X or self.type == EdgeType.DIRECTED_TO_Y

    def get_source(self):
        if self.type == EdgeType.DIRECTED_TO_X:
            return self.y
        if self.type == EdgeType.DIRECTED_TO_Y:
            return self.x
        raise RuntimeError("Edge is not directed")

    def get_target(self):
        if self.type == EdgeType.DIRECTED_TO_X:
            return self.x
        if self.type == EdgeType.DIRECTED_TO_Y:
            return self.y
        raise RuntimeError("Edge is not directed")

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.type == other.type

    def __hash__(self):
        return hash(self.x) ^ (hash(self.y) << 1) ^ (hash(self.type) << 2)


class EdgeQueueSet:
    def __init__(self):
        self.edges_queue = deque()
        self.edges_set = set()

    def push_directed(self, x, y):
        edge = Edge(x, y, EdgeType.DIRECTED_TO_Y)
        if edge not in self.edges_set:
            self.edges_queue.append(edge)
            self.edges_set.add(edge)

    def push_undirected(self, x, y):
        if x > y:
            x, y = y, x
        edge = Edge(x, y, EdgeType.UNDIRECTED)
        if edge not in self.edges_set:
            self.edges_queue.append(edge)
            self.edges_set.add(edge)

    def pop(self):
        edge = self.edges_queue.popleft()
        self.edges_set.remove(edge)
        return edge

    def empty(self):
        return len(self.edges_queue) == 0


class EdgeModification:
    def __init__(self, x, y, old_type, new_type):
        self.x = x
        self.y = y
        self.old_type = old_type
        self.new_type = new_type

    def is_reverse(self):
        return (
            self.old_type == EdgeType.DIRECTED_TO_X and self.new_type == EdgeType.DIRECTED_TO_Y
        ) or (self.old_type == EdgeType.DIRECTED_TO_Y and self.new_type == EdgeType.DIRECTED_TO_X)

    def is_new_directed(self):
        return self.new_type == EdgeType.DIRECTED_TO_X or self.new_type == EdgeType.DIRECTED_TO_Y

    def is_old_directed(self):
        return self.old_type == EdgeType.DIRECTED_TO_X or self.old_type == EdgeType.DIRECTED_TO_Y

    def is_new_undirected(self):
        return self.new_type == EdgeType.UNDIRECTED

    def is_old_undirected(self):
        return self.old_type == EdgeType.UNDIRECTED

    def get_new_target(self):
        if self.new_type == EdgeType.DIRECTED_TO_X:
            return self.x
        if self.new_type == EdgeType.DIRECTED_TO_Y:
            return self.y
        raise RuntimeError("Edge is not directed")

    def get_new_source(self):
        if self.new_type == EdgeType.DIRECTED_TO_X:
            return self.y
        if self.new_type == EdgeType.DIRECTED_TO_Y:
            return self.x
        raise RuntimeError("Edge is not directed")

    def get_old_target(self):
        if self.old_type == EdgeType.DIRECTED_TO_X:
            return self.x
        if self.old_type == EdgeType.DIRECTED_TO_Y:
            return self.y
        raise RuntimeError("Edge is not directed")

    def get_old_source(self):
        if self.old_type == EdgeType.DIRECTED_TO_X:
            return self.y
        if self.old_type == EdgeType.DIRECTED_TO_Y:
            return self.x
        raise RuntimeError("Edge is not directed")

    def get_modification_id(self):
        if self.old_type == EdgeType.NONE and self.is_new_undirected():
            return 1
        if self.old_type == EdgeType.NONE and self.is_new_directed():
            return 2
        if self.is_old_undirected() and self.new_type == EdgeType.NONE:
            return 3
        if self.is_old_undirected() and self.is_new_directed():
            return 4
        if self.is_old_directed() and self.new_type == EdgeType.NONE:
            return 5
        if self.is_old_directed() and self.is_new_undirected():
            return 6
        if self.is_reverse():
            return 7
        raise RuntimeError("Invalid modification")

    def __repr__(self):
        return (
            f"{self.x} {EdgeType.to_string(self.old_type)} {self.y} "
            f"becomes {self.x} {EdgeType.to_string(self.new_type)} {self.y}"
        )


class EdgeModificationsMap:
    def __init__(self):
        self.edge_modifications = {}

    def update_edge_directed(self, x, y, old_type):
        if x > y:
            x, y = y, x
            new_type = EdgeType.DIRECTED_TO_X
            if old_type == EdgeType.DIRECTED_TO_X:
                old_type = EdgeType.DIRECTED_TO_Y
        else:
            new_type = EdgeType.DIRECTED_TO_Y
        self.update_edge_modification(x, y, old_type, new_type)

    def update_edge_undirected(self, x, y, old_type):
        if x > y:
            x, y = y, x
            if old_type == EdgeType.DIRECTED_TO_X:
                old_type = EdgeType.DIRECTED_TO_Y
            elif old_type == EdgeType.DIRECTED_TO_Y:
                old_type = EdgeType.DIRECTED_TO_X
        self.update_edge_modification(x, y, old_type, EdgeType.UNDIRECTED)

    def update_edge_none(self, x, y, old_type):
        if x > y:
            x, y = y, x
            if old_type == EdgeType.DIRECTED_TO_X:
                old_type = EdgeType.DIRECTED_TO_Y
            elif old_type == EdgeType.DIRECTED_TO_Y:
                old_type = EdgeType.DIRECTED_TO_X
        self.update_edge_modification(x, y, old_type, EdgeType.NONE)

    def update_edge_modification(self, small, big, old_type, new_type):
        edge_key = (small, big)
        if edge_key in self.edge_modifications:
            oldest_type = self.edge_modifications[edge_key].old_type
            if oldest_type == new_type:
                del self.edge_modifications[edge_key]
            else:
                self.edge_modifications[edge_key].new_type = new_type
        else:
            self.edge_modifications[edge_key] = EdgeModification(small, big, old_type, new_type)

    def clear(self):
        self.edge_modifications.clear()

    def __iter__(self):
        return iter(self.edge_modifications.values())
