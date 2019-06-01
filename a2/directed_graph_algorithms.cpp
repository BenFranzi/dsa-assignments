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

////iterator pattern
//for (auto item : d) {
//std::cout << item << " --> ";
//for(auto neighbour = d.nbegin(item); neighbour != d.nend(item); ++neighbour) {
//std::cout << *neighbour;
//}
//std::cout << std::endl;
//}



//    pending.push(v);
//    pendingArray.push_back(v);
//    for (auto neighbour = d.nbegin(v); neighbour != d.nend(v); ++neighbour) {
//        for (vertex item : pendingArray) {
//            if (item == *neighbour) //duplicate item in pending;
//                return false;
//            DFS(d, *neighbour, pending, pendingArray);
//        }
//    }
//    auto popped = pending.pop();
//    pendingArray.erase(popped);
//    return true;


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

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
    return std::list<vertex>();
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
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

/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
    return std::unordered_map<vertex, std::size_t>();
}

