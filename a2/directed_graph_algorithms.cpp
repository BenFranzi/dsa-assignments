#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

/**
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one
 * vertex, and returns to it.
 */
template <typename vertex> bool is_dag(const directed_graph<vertex> & d) {
    std::unordered_set<vertex> visited;
    std::unordered_set<vertex> pending;
    //Iterate over each vertex in the graph
    for (auto vert : d) {
        //Clear the visit and pending sets for next iteration
        visited.clear();
        pending.clear();
        //If recursive call returns having found a cycle it is not acyclic
        if (is_cycle(d, vert, visited, pending))
            return false;
    }
    //If no cycle is found then the graph is acyclic
    return true;
}

/**
 * Determines if a cycle is found under a given vertex
 * @param d - reference to the whole graph
 * @param v - current vertex to iterate from
 * @param visited - list of the vertices that have been visited
 * @param pending - list of vertices visited on current path
 * @return true if a cycle is found, false if the no cycle is found
 */
template <typename vertex> bool is_cycle(const directed_graph<vertex> &d, const vertex &v, std::unordered_set<vertex> & visited, std::unordered_set<vertex> & pending) {
    //If the cycle has not been visited
    if (visited.count(v) == 0) {
        visited.insert(v);
        pending.insert(v);
        //For each neighbour of the given vertex
        for(auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
            auto curr_neigh = *neighbour;
            // If the current neighbour has not be visited and the recursive cycle has found a cycle then a cycle exists
            if (visited.count(*neighbour) == 0 && is_cycle(d, *neighbour, visited, pending))
                return true;
            //If the current neighbour already exists in the pending queue then a cycle exists
            else if (pending.count(*neighbour) != 0)
                return true;
        }
    }
    //Remove the current vertex from the pending as no cycle was found from that point
    pending.erase(v);
    return false;
}

template <typename vertex> void sort(const directed_graph<vertex> &d, const vertex &v, std::unordered_set<vertex> & visited, std::stack<vertex> & sorted) {
    visited.insert(v);
    for(auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
        bool exists_in_visited = visited.count(*neighbour) != 0;
        if (!exists_in_visited)
            sort(d, *neighbour, visited, sorted);
    }
    sorted.push(v);
}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
    std::unordered_set<vertex> visited;
    std::stack<vertex> sorted;

    for (auto vert : d) {
        bool exists_in_visited = visited.count(vert) != 0;
        if (!exists_in_visited) {
            sort(d, vert, visited, sorted);
        }
    }

    std::list<vertex> result;
    while (!sorted.empty())
    {
        result.push_back(sorted.top());
        sorted.pop();
    }

    return result;
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
    //Check if graph is empty or only contains a single element
    if (d.num_vertices() < 1) {
        return true;
    }

    std::list<vertex> sorted = topological_sort(d);
    std::vector<vertex> topology;

    for (auto vert : sorted) {
        topology.push_back(vert);
    }

    for (int i=0; i<topology.size() - 1; ++i) {
        if (!(d.adjacent(topology[i], topology[i+1]))) {
            return false;
        }
    }
    return true;
}

template <typename vertex> void convert_graph(const directed_graph<vertex> & d, directed_graph<vertex> & graph, vertex & v, std::unordered_set<vertex> & visited) {
    if (visited.count(v) == 0) {
        visited.insert(v);
        for(auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
            graph.add_edge(v, *neighbour);
            graph.add_edge(*neighbour, v);
        }
    }
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
    directed_graph<vertex> graph;
    std::unordered_set<vertex> visited;
    for (auto vert : d) {
        graph.add_vertex(vert);
    }
    for (auto vert : d) {
        convert_graph(d, graph, vert, visited);
    }
    return strongly_connected_components(graph);
}

template <typename vertex> void tarjan(
        const directed_graph<vertex> & d,
        const vertex & v,
        std::unordered_map<vertex, int> & index,
        std::unordered_map<vertex, int> & low,
        std::stack<vertex> & st,
        std::unordered_set<vertex> & onstack,
        std::vector<std::vector<vertex>> & result
        ) {
    static int time = 0;

    index[v] = time;
    low[v] = time;
    time++;
    st.push(v);
    onstack.insert(v);

    for(auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
        if (index.count(*neighbour) == 0) {
            tarjan(d, *neighbour, index, low, st, onstack, result);
            low[v] = (low[v] < low[*neighbour]) ? low[v] : low[*neighbour];
        } else if (onstack.count(*neighbour) != 0) {
            low[v] = (low[v] < index[*neighbour]) ? low[v] : index[*neighbour];
        }
    }

    if (low[v] == index[v]) {
        std::vector<vertex> new_components;
        auto w = st.top();
        do {
            w = st.top();
            onstack.erase(w);
            new_components.push_back(w);
            st.pop();
        } while (w != v);
        result.push_back(new_components);
    }
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */
template <typename vertex> std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
    std::unordered_map<vertex, int> index;
    std::unordered_map<vertex, int> low;
    std::stack<vertex> st;
    std::unordered_set<vertex> onstack;
    std::vector<std::vector<vertex>> result;

    for (auto vert : d) {
        if (index.count(vert) == 0) {
            tarjan(d, vert, index, low, st, onstack, result);
        }
    }

    return result;
}


template <typename vertex> void find_distance(const directed_graph<vertex> & d, const vertex & u, std::unordered_map<vertex, int> & results, const int & count, std::unordered_set<vertex> & visited) {

    for(auto neighbour = d.nbegin(u); neighbour != d.nend(u); ++neighbour) {
        auto curr = *neighbour;
        int curr_count = results[*neighbour];
        if  (count < curr_count) {
            results[*neighbour] = count;
            find_distance(d, *neighbour, results, count + 1, visited);
        }
    }
}
/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, int> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
    std::unordered_map<vertex, int> results;
    std::unordered_set<vertex> visited;
//    int count = 0;
    for (auto vert : d) {
        if (vert == u) {
            results[vert] = 0;
        } else {
            results[vert] = d.num_vertices() + 1;
        }
    }

    //static_cast<int>
    find_distance(d, u, results, 1, visited);
    return results;
}

