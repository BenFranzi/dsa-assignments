cmake_minimum_required(VERSION 3.14)
project(a2)

set(CMAKE_CXX_STANDARD 14)

add_executable(a2 main.cpp directed_graph.hpp directed_graph_algorithms.cpp)