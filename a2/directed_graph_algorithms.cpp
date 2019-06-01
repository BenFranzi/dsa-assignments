/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
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


//Print relationships
template <typename vertex> void print_graph(const directed_graph<vertex> & d) {
    for (auto item : d) {
        std::cout << item << " --> ";
        for(auto neighbour = d.nbegin(item); neighbour != d.nend(item); ++neighbour) {
            std::cout << *neighbour;
        }
        std::cout << std::endl;
    }
}

//DFS
template <typename vertex> bool is_cycle(const directed_graph<vertex> &d, const vertex &v, std::unordered_set<vertex> & visited, std::unordered_set<vertex> & pending) {
    if (visited.count(v) == 0) {
        visited.insert(v);
        pending.insert(v);
        for(auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
            auto curr_neigh = *neighbour;
            bool exists_in_visited = visited.count(*neighbour) != 0;
            bool exists_in_pending = pending.count(*neighbour) != 0;
            if (!exists_in_visited && is_cycle(d, *neighbour, visited, pending))
                return true;
            else if (exists_in_pending)
                return true;
        }
    }
    pending.erase(v);
    return false;
}
//END HELPERS

/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one
 * vertex, and returns to it.
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
    std::unordered_set<vertex> visited;
    std::unordered_set<vertex> pending;
    for (auto vert : d) {
        visited.clear();
        pending.clear();
        bool x = is_cycle(d, vert, visited, pending);
        if (x)
            return false;
    }
    return true;
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

template <typename vertex> bool is_path(const directed_graph<vertex> &d, const std::size_t size, const vertex &v) {
//    if (visited.count(v) == 0) {
//        visited.insert(v);
//        pending.insert(v);
//        for(auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
//            auto curr_neigh = *neighbour;
//            bool exists_in_visited = visited.count(*neighbour) != 0;
//            bool exists_in_pending = pending.count(*neighbour) != 0;
//            if (!exists_in_visited && is_cycle(d, *neighbour, visited, pending))
//                return true;
//            else if (exists_in_pending)
//                return true;
//        }
//    }
//    pending.erase(v);
    return false;

}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
    std::unordered_set<vertex> visited;

    for (auto vert : d) {
        bool x = is_path(d, d.num_vertices(), vert);
        if (x)
            return true;
    }
    return false;
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
    return std::vector<std::vector<vertex>>();
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
    return std::vector<std::vector<vertex>>();
}


template <typename vertex> void find_distance(const directed_graph<vertex> & d, const vertex & u, std::unordered_map<vertex, int> & results, const int & count, std::unordered_set<vertex> & visited) {

    for(auto neighbour = d.nbegin(u); neighbour != d.nend(u); ++neighbour) {
        auto curr = *neighbour;
        int curr_count = results[*neighbour];
        if (curr_count == -1) {
            results[*neighbour] = count;
            find_distance(d, *neighbour, results, count + 1, visited);
        } else if (curr_count > count) {
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

