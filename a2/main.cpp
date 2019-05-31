#include <iostream>
#include "directed_graph_algorithms.cpp"
#include "directed_graph.hpp"

int main() {
    std::cout << "Hello, World!" << std::endl;
    directed_graph<char> graph;

    char a = 'a';
    char b = 'b';
    char c = 'c';
    char d = 'd';
    char e = 'e';
    char f = 'f';
    char g = 'g';
    char h = 'h';


    graph.add_vertex(a);
    graph.add_vertex(b);
    graph.add_vertex(c);

    graph.add_edge(a, b);
    graph.add_edge(a, c);

    directed_graph<char> empty;

    std::cout << is_dag(graph) << std::endl;

    return 0;
}