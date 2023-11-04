#include <iostream>
#include "PDAG.h" // Include the header file for the PDAG class

int main() {
    // Create a PDAG instance with nodes 0, 1, 2, 3, and 4
    PDAG pdag(5);

    pdag.add_directed_edge(0, 1);
    pdag.add_undirected_edge(2, 3);

    std::cout << "Number of edges: " << pdag.get_number_of_edges() << std::endl;
    // show the list of directed edges
    std::cout << "Directed edges: ";
    for (const auto &edge: pdag.get_directed_edges()) {
        std::cout << "(" << edge.first << ", " << edge.second << ") ";
    }

    std::cout << std::endl << "Is clique (1,2,3): ";
    std::cout << pdag.is_clique({1, 2, 3}) << std::endl;

    return 0;
}