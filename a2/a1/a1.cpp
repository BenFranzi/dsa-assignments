#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
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

template <typename vertex> class vertex_iterator;
template <typename vertex> class neighbour_iterator;
template <typename vertex> class directed_graph;


template <typename vertex>
class directed_graph {

private:
    std::size_t edges_count; //Holds the number of edges in the graph
    int get_position(const vertex&) const; //Gets the index of a vertex from its pointer
    int is_valid(int) const; //Checks if get_position returned a result
    std::vector<vertex> vertices; //Holds the vertices in the graph
    std::vector<std::vector<bool>> adj_matrix; //Adjacency matrix that will hold either a True or a False for a directed arc

public:
    directed_graph(); //A constructor for directed_graph. The graph should start empty.
    ~directed_graph(); //A destructor. Depending on how you do things, this may
    //not be necessary.

    bool contains(const vertex&) const; //Returns true if the given vertex is in the graph, false otherwise.

    bool adjacent(const vertex&, const vertex&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

    void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
    void add_edge(const vertex&, const vertex&); //Adds an edge from the first vertex to the second.

    void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
    void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.

    std::size_t in_degree(const vertex&) const; //Returns number of edges coming in to a vertex.
    std::size_t out_degree(const vertex&) const; //Returns the number of edges leaving a vertex.
    std::size_t degree(const vertex&) const; //Returns the degree of the vertex (both in and out edges).

    std::size_t num_vertices() const; //Returns the total number of vertices in the graph.
    std::size_t num_edges() const; //Returns the total number of edges in the graph.

    std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
    std::vector<vertex> const_get_vertices() const; //Returns a const vector containing all the vertices for the vertex_iterator.
    std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

    vertex_iterator<vertex> begin(); //Returns a vertex_iterator pointing to the start of the vertex set.

    vertex_iterator<vertex> end(); //Returns a vertex_iterator pointing to one-past-the-end of the vertex set.

    neighbour_iterator<vertex> nbegin(const vertex&); //Returns a neighbour_iterator pointing to the start of the neighbour set for the given vertex.
    neighbour_iterator<vertex> nend(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end of the neighbour set for the given vertex.

    std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
    std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

    directed_graph<vertex> out_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges.
    directed_graph<vertex> in_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the in-edges.

    bool reachable(const vertex&, const vertex&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
};

//The vertex_iterator class provides an iterator
//over the vertices of the graph.
//This is one of the harder parts, so if you're
//not too comfortable with C++ leave this for last.
//If you are, there are many ways of doing this,
//as long as it passes the tests, it's okay.
//You may want to watch the videos on iterators before starting.
template <typename vertex>
class vertex_iterator {

private:
    const std::vector<vertex> vertices; //Graphs vertices
    std::size_t position; //Current vertex position in the vertices vector

public:
    vertex_iterator(const vertex_iterator<vertex>&);
    vertex_iterator(const directed_graph<vertex>&, std::size_t);
    ~vertex_iterator();
    vertex_iterator<vertex> operator=(const vertex_iterator<vertex>&);
    bool operator==(const vertex_iterator<vertex>&) const;
    bool operator!=(const vertex_iterator<vertex>&) const;
    vertex_iterator<vertex> operator++();
    vertex_iterator<vertex> operator++(int);
    vertex operator*();
    vertex* operator->();
};

//The neighbour_iterator class provides an iterator
//over the neighbours of a given vertex. This is
//probably harder (conceptually) than the graph_iterator.
//Unless you know how iterators work.
template <typename vertex>
class neighbour_iterator {

private:
    directed_graph<vertex> graph; //The graph to be iterated
    const vertex curr; //Current vertex
    std::size_t position; //Current neighbour vertex index

public:
    neighbour_iterator(const neighbour_iterator<vertex>&);
    neighbour_iterator(const directed_graph<vertex>&, const vertex&, std::size_t);
    ~neighbour_iterator();
    neighbour_iterator<vertex> operator=(const neighbour_iterator<vertex>&);
    bool operator==(const neighbour_iterator<vertex>&) const;
    bool operator!=(const neighbour_iterator<vertex>&) const;
    neighbour_iterator<vertex> operator++();
    neighbour_iterator<vertex> operator++(int);
    vertex operator*();
    vertex* operator->();
};


//Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

/**
 * Directed graph constructor
 * Default vertices and edge counts to zero
 */
template <typename vertex> directed_graph<vertex>::directed_graph() {
    edges_count = 0;
}

/**
 * Destructor
 */
template <typename vertex> directed_graph<vertex>::~directed_graph() {}

/**
 * Checks if graph contains vertex u
 * @param u
 * @return bool - if vertex exists in graph
 */
template <typename vertex> bool directed_graph<vertex>::contains(const vertex& u) const {
    //Iterate over vertices
    for (vertex v : vertices) {
        if (u == v) {
            //Return true if found
            return true;
        }
    }
    //Return false if not found
    return false;
}

/**
 * Checks if v is an out neighbour of u
 * @param u - start vertex
 * @param v - end vertex
 * @return if vertex v is out edge of u
 */
template <typename vertex> bool directed_graph<vertex>::adjacent(const vertex& u, const vertex& v) const {
    int u_pos = get_position(u);
    int v_pos = get_position(v);
    //If two valid positions are given
    if (is_valid(u_pos) && is_valid(v_pos)) {
        //Return value from adjacency matrix
        return adj_matrix[u_pos][v_pos];
    }
    return false;
}

/**
 * Adds a vertex to the graph
 * @param u - pointer to vertex being added
 */
template <typename vertex> void directed_graph<vertex>::add_vertex(const vertex& u) {
    //Check if the vertex already exists
    if (!is_valid(get_position(u))) {
        //Add new vertex to the list
        vertices.push_back(u);
        //Add new vertex to end of rows in matrix
        for (auto& row : adj_matrix) {
            row.push_back(false);
        }
        //Add new row of false at the bottom of the matrix
        adj_matrix.push_back(std::vector<bool>(vertices.size(), false));
    }
}

/**
 * Add edge between two vertices in graph
 * @param u - start vertex
 * @param v - end vertex
 */
template <typename vertex> void directed_graph<vertex>::add_edge(const vertex& u, const vertex& v) {
    // Add arc from u to v
    int u_pos = get_position(u);
    int v_pos = get_position(v);

    // If valid vertices and there is no existing relationship
    if (is_valid(u_pos) && is_valid(v_pos) && !adj_matrix[u_pos][v_pos]) {
        //Set relationship in adjacency matrix
        adj_matrix[u_pos][v_pos] = true;
        ++edges_count;
    }
}

/**
 * Remove vertex from graph and all corresponding relationships
 * @param u - vertex to be removed
 */
template <typename vertex> void directed_graph<vertex>::remove_vertex(const vertex& u) {
    int pos = get_position(u);
    //Check if vertex exists in graph
    if (is_valid(pos)) {
        // Remove row for vertices
        adj_matrix.erase(adj_matrix.begin() + pos);
        // Remove column of vertices
        for (auto& row : adj_matrix) {
            row.erase(row.begin() + pos);
        }
        // Remove vertex from list
        vertices.erase(vertices.begin() + pos);
    }
}

/**
 * Removes edge between two vertices
 * @param u - start vertex
 * @param v - end vertex
 */
template <typename vertex> void directed_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
    int u_pos = get_position(u);
    int v_pos = get_position(v);
    // If valid vertices and there is an existing relationship
    if (is_valid(u_pos) && is_valid(v_pos) && adj_matrix[u_pos][v_pos]) {
        //Set relationship in adjacency matrix
        adj_matrix[u_pos][v_pos] = false;
        --edges_count;
    }
}

/**
 * Returns the in degree of a given vertex
 * @param u - vertex
 * @return degree of in edges of vertex
 */
template <typename vertex> std::size_t directed_graph<vertex>::in_degree(const vertex& u) const {
    int count = 0;
    int pos = get_position(u);
    if (is_valid(pos)) {
        //Iterate over "in" edges for vertex
        //Iterates down the column in the adjacency matrix
        for (std::vector<bool> row : adj_matrix) {
            //If true in column for adjacency matrix, increment count
            if (row[pos]) ++count;
        }
    }
    return count;
}

/**
 * Returns the out degree of a given vertex
 * @param u - vertex
 * @return degree of out edges of vertex
 */
template <typename vertex> std::size_t directed_graph<vertex>::out_degree(const vertex& u) const {
    int count = 0;
    int pos = get_position(u);
    if (is_valid(pos)) {
        //Iterate over "out" edges for vertex
        //Iterates across the row in the adjacency matrix
        for (bool out : adj_matrix[pos]) {
            //If true in row for adjacency matrix, increment count
            if (out) ++count;
        }
    }
    return count;
}

/**
 * Get total degree of a given vertex
 * @param u - vertex
 * @return total in and out degree of vertex
 */
template <typename vertex> std::size_t directed_graph<vertex>::degree(const vertex& u) const {
    return in_degree(u) + out_degree(u);
}

/**
 * Return total vertices in graph;
 * @tparam vertex
 * @return
 */
template <typename vertex> std::size_t directed_graph<vertex>::num_vertices() const { return vertices.size(); }


/**
 * Return the total number of edges in the graph
 * @return total edges in graph
 */
template <typename vertex> std::size_t directed_graph<vertex>::num_edges() const { return edges_count; }

/**
 * Returns all vertices in graph
 * @return vertices in graph
 */
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_vertices() { return vertices; }


/**
 * Returns all vertices in graph as const for vertex_iterator
 * @return vertices in graph as const
 */
template <typename vertex> std::vector<vertex> directed_graph<vertex>::const_get_vertices() const { return vertices; }


/**
 * Returns out neighbours at a given vertex
 * @param u - vertex
 * @return out neighbours of u
 */
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_neighbours(const vertex& u) {
    int pos = get_position(u);
    std::vector<vertex> neighbours;
    if (is_valid(pos)) {
        //Check adjacency matrix for all vertices
        for (int i=0; i<num_vertices(); ++i) {
            //If out edge found, add vertex to out vector
            //Checks for out edges across row
            if (adj_matrix[pos][i]) {
                neighbours.push_back(vertices[i]);
            }
        }
    }
    return neighbours;
}

/**
 * Begin
 * @return returns vertex iterator at start of vertices array
 */
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::begin() {
    return vertex_iterator(*this, 0);
}

/**
 * End
 * @return returns vertex iterator at end of vertices array
 */
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::end() {
    return vertex_iterator(*this, vertices.size());
}

/**
 * Begin neighbours
 * @param u - vertex to iterate
 * @return returns iterator at start for a given vertex's neighbour vertices
 */
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nbegin(const vertex& u) {
    return neighbour_iterator(*this, u, 0);
}

/**
 * End neighbours
 * @param u - vertex to iterate
 * @return returns iterator at end for a given vertex's neighbour vertices
 */
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nend(const vertex& u) {
    return neighbour_iterator(*this, u, this->get_neighbours(u).size());
}

/**
 * Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex
 * @param u - start vertex
 * @return vector of vertices that were found after the depth first traversal at that given vertex
 */
template <typename vertex> std::vector<vertex> directed_graph<vertex>::depth_first(const vertex& u) {
    std::vector<bool> visited(this->num_vertices(), false);
    //Using a stack for depth first traversal
    std::stack<vertex> pending;
    std::vector<vertex> result;

    //Push start to pending list
    pending.push(u);

    //Repeat until all pending items have been handled
    while (!pending.empty()) {
        //Get the current vertex and position in array
        vertex curr = pending.top();
        int curr_pos = this->get_position(curr);
        //Remove current from pending
        pending.pop();
        //If the vertex hasn't been visited yet
        if (!visited[curr_pos]) {
            visited[curr_pos] = true;
            //Add vertex to list of visited vertices
            result.push_back(curr);
            //Iterate over out edges of current vertex
            for (int i=this->num_vertices(); i > 0; --i) {
                if (adj_matrix[curr_pos][i-1]) {
                    //If out edge, add to pending if edge exists
                    pending.push(vertices[i-1]);
                }
            }
        }
    }
    return result;
}

/**
 * Returns the vertices of the graph in the order they are visited in by a breadth-first traversal starting at the given vertex
 * @param u - start vector
 * @return vector of vertices that were found after the breadth first traversal at that given vertex
 */
template <typename vertex> std::vector<vertex> directed_graph<vertex>::breadth_first(const vertex& u) {
    std::vector<bool> visited(this->num_vertices(), false);
    //Using a queue for breadth first traversal
    std::queue<vertex> pending;
    std::vector<vertex> result;

    //Push start to pending list
    pending.push(u);

    //Repeat until all pending items have been handled
    while(!pending.empty()) {
        //Get the current vertex and position in array
        vertex curr = pending.front();
        int curr_pos = this->get_position(curr);
        //Remove current from pending
        pending.pop();
        if (!visited[curr_pos]) {
            visited[curr_pos] = true;
            //Add vertex to list of visited vertices
            result.push_back(curr);
            //Iterate over out edges of current vertex
            for (int i=0; i<this->num_vertices(); ++i) {
                if (adj_matrix[curr_pos][i]) {
                    //If out edge, add to pending if edge exists
                    pending.push(vertices[i]);
                }
            }
        }
    }
    return result;
}

/**
 * Returns a spanning tree of the graph starting at the given vertex using the out-edges
 * @param u - start vertex
 * @return spanning tree of out edges from vertex
 */
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::out_tree(const vertex& u) {
    std::vector<bool> visited(this->num_vertices(), false);
    //Implemented spanning tree using depth first traversal so a stack is used
    std::stack<vertex> pending;
    directed_graph<vertex> tree;

    //Add start vertex
    pending.push(u);
    tree.add_vertex(u);

    //Repeat until all pending items have been handled
    while(!pending.empty()) {
        //Get the current vertex and position in array
        vertex curr = pending.top();
        int curr_pos = get_position(curr);
        //Remove current from pending
        pending.pop();
        if (!visited[curr_pos]) {
            //Add vertex to list of visited vertices
            visited[curr_pos] = true;
            //Iterate over out edges of current vertex
            for (int i = 0; i < num_vertices(); ++i) {
                //If out edge exists and vertex hasn't been made already
                if (adj_matrix[curr_pos][i] && !tree.contains(vertices[i])) {
                    //Add the vertex and make the edge relationship
                    tree.add_vertex(vertices[i]);
                    tree.add_edge(curr, (vertices[i]));

                    //Add the found neighbour to the pending list
                    pending.push(vertices[i]);
                }
            }
        }
    }
    return tree;
}

/**
 * Returns a spanning tree of the graph starting at the given vertex using the in-edges
 * @param u - start vertex
 * @return spanning tree of in edges from vertex
 */
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::in_tree(const vertex& u) {
    std::vector<bool> visited(this->num_vertices(), false);
    //Implemented spanning tree using depth first traversal so a stack is used
    std::stack<vertex> pending;
    directed_graph<vertex> tree;

    //Add start vertex
    pending.push(u);
    tree.add_vertex(u);

    //Repeat until all pending items have been handled
    while(!pending.empty()) {
        //Get the current vertex and position in array
        vertex curr = pending.top();
        int curr_pos = get_position(curr);
        //Remove current from pending
        pending.pop();
        if (!visited[curr_pos]) {
            //Add vertex to list of visited vertices
            visited[curr_pos] = true;
            //Iterate over out edges of current vertex
            for (int i=0; i<num_vertices(); ++i) {
                //If in edge exists and vertex hasn't been made already
                if (adj_matrix[i][curr_pos] && !tree.contains(vertices[i])) {
                    //Add the vertex and make the edge relationship backwards
                    tree.add_vertex(vertices[i]);
                    tree.add_edge(vertices[i], curr);

                    //Add the found neighbour to the pending list
                    pending.push(vertices[i]);
                }
            }
        }
    }
    return tree;
}

/**
 * Determines if v is reachable from v
 * @param u - start vertex
 * @param v - end vertex
 * @return returns true if v is reachable from u
 */
template <typename vertex> bool directed_graph<vertex>::reachable(const vertex& u, const vertex& v) const {
    directed_graph<vertex> curr = *this;
    //Create a spanning tree from u vertex
    directed_graph<vertex> tree = curr.out_tree(u);
    //Return true if the spanning tree added vertex v to the tree
    return tree.contains(v);
}


/***************************
  VERTEX ITERATOR
****************************/

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const vertex_iterator<vertex>& other) : vertices(other.vertices), position(other.position) {}
template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const directed_graph<vertex>& graph, std::size_t position) : vertices(graph.const_get_vertices()), position(position) {}
template <typename vertex> vertex_iterator<vertex>::~vertex_iterator() {}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator=(const vertex_iterator<vertex>& other) {
    vertices = other.vertices;
    position = other.position;
}

template <typename vertex> bool vertex_iterator<vertex>::operator==(const vertex_iterator<vertex>& other) const {
    return (vertices == other.vertices && position == other.position);
}

template <typename vertex> bool vertex_iterator<vertex>::operator!=(const vertex_iterator<vertex>& other) const {
    return (vertices != other.vertices || position != other.position);
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++() {
    ++position;
    return *this;
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++(int) {
    auto temp = this;
    ++position;
    return temp;
}

template <typename vertex> vertex vertex_iterator<vertex>::operator*() {
    return vertices[position];
}

template <typename vertex> vertex* vertex_iterator<vertex>::operator->() {
    return vertices + position;
}


/***************************
  NEIGHBOUR ITERATOR
****************************/

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const neighbour_iterator<vertex>& other) : graph(other.graph), curr(other.curr), position(other.position) {}
template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const directed_graph<vertex>& graph, const vertex& u, std::size_t position) : graph(graph), curr(u), position(position) {}
template <typename vertex> neighbour_iterator<vertex>::~neighbour_iterator() {}
template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator=(const neighbour_iterator<vertex>& other) {
    graph = other.graph;
    curr = other.curr;
    position = other.position;
}
template <typename vertex> bool neighbour_iterator<vertex>::operator==(const neighbour_iterator<vertex>& other) const {
    return (curr == other.curr && position == other.position);
}
template <typename vertex> bool neighbour_iterator<vertex>::operator!=(const neighbour_iterator<vertex>& other) const {
    return (curr != other.curr || position != other.position);
}
template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++() {
    ++position;
    return *this;
}
template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++(int) {
    auto temp = this;
    ++position;
    return temp;
}
template <typename vertex> vertex neighbour_iterator<vertex>::operator*() {
    return graph.get_neighbours(curr)[position];
}
template <typename vertex> vertex* neighbour_iterator<vertex>::operator->() {
    return &graph.get_neighbours(curr)[position];
}


/***************************
  PRIVATE
****************************/

/**
 * Get the position of vertex in vector from the pointer
 * @param v - vertex
 * @return 0 to n if found, -1 if not found
 */
template <typename vertex> int directed_graph<vertex>::get_position(const vertex& v) const {
    //TODO: see if there is a better method of determining the index
    for (int i=0; i < vertices.size(); i++) {
        if (vertices[i] == v) {
            return i;
        }
    }
    // Failed to find in array
    return -1;
}

/**
 * Checks if the index is valid
 * @param pos - index value
 * @return if the index is valid
 */
template <typename vertex> int directed_graph<vertex>::is_valid(int pos) const {
    //Checks if the get_position call returned a value that was valid
    return pos >= 0;
    //This function exists to make the code more readable
}

#endif