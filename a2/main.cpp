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
    graph.add_vertex(d);
    graph.add_vertex(e);
    graph.add_vertex(f);

    graph.add_edge(a, b);
    graph.add_edge(b, e);
    graph.add_edge(a, c);
    graph.add_edge(c, d);
    graph.add_edge(d, e);
    graph.add_edge(d, b);
    graph.add_edge(c, f);

    ////CYCLE
//    graph.add_edge(b, c);
    /////////
    print_graph(graph);
    std::cout << is_dag(graph, a) << std::endl;

    return 0;
}