#Parallelization for Locality-Based Algos
##Graph Algorithms
Here, graph algortihms are those locality-based overlapping community detection algorithms.

##Source Code
- [src/algorithm](src/algorithm), implemented graph algorithms
- [src/algorithm_demo](src/algorithm_demo), simple demo source codes to execute graph algorithm program
- [src/evaluation_metrics](src/evaluation_metrics), evaluation metrics for graph algorithms
- [src/parallel_utils](src/parallel_utils), parallel utilities for accelerating computations of graph algorithms
- [src/util](src/util), utilities for graph input and pretty printing

##Demo Files
- [demo_files/demo_graph.csv](demo_files/demo_graph.csv), toy graph in edge-list format
- [demo_files/demo_result.txt](demo_files/demo_result.txt),
computation result of connected-iterative-scan algorithm, where each line is a community

##Dataset
- [small_datasets/collaboration_edges_input.csv](small_datasets/collaboration_edges_input.csv), from uci repo
- [small_datasets/karate_edges_input.csv](small_datasets/karate_edges_input.csv), from uci repo
