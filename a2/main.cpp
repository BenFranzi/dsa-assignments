#include <iostream>
#include "directed_graph_algorithms.cpp"
#include "directed_graph.hpp"

int main() {
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
//    graph.add_vertex(f);

//    graph.add_edge(a, b);
//    graph.add_edge(b, e);
//    graph.add_edge(a, c);
//    graph.add_edge(c, d);
//    graph.add_edge(d, e);
//    graph.add_edge(d, b);
//    graph.add_edge(c, f);

    //TEST DAG
    //print_graph(graph);
    //std::cout << "should be true  (1): " << is_dag(graph) << std::endl;
    //graph.add_edge(b, c);
    //std::cout << "should be false (0): " << is_dag(graph) << std::endl;

    //TEST TOPO
    //std::cout << "true(1) || false(0) " << is_hamiltonian_dag(graph) << std::endl;

    //TEST SHORTEST
    //graph.add_edge(a, b);
    //graph.add_edge(b, d);
    //graph.add_edge(b, e);
    //graph.add_edge(a, e);
    //graph.add_edge(d, c);
    //shortest_distances(graph, b);

    //TEST HAMILTONIAN
    graph.add_edge(a, b);
    graph.add_edge(a, c);

    graph.add_edge(b, e);
    graph.add_edge(b, d);
    graph.add_edge(c, d);
    graph.add_edge(d, e);

    std::cout << "should be false (0): " << is_hamiltonian_dag(graph) << std::endl;
    graph.add_edge(c, b);
    std::cout << "should be true  (1): " << is_hamiltonian_dag(graph) << std::endl;
    return 0;
}