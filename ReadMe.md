# Parallelization for Locality-Based Algos

## Issues
In community merging part for cis, best strategy is not applied, since cis expect totally including relationship.
Extra cost is introduced in current implementation, a better solution can be that, build index for large community,
iterate through small community and judge whether v in small community is found in the large community.

## Graph Algorithms
Here, graph algorithms are those locality-based overlapping community detection algorithms.

## Source Code

content | detail
--- | ---
[src/algorithm](src/algorithm) | implemented graph algorithms
[src/algorithm_demo](src/algorithm_demo) | simple demo source codes to execute graph algorithm program
[src/parallel_utils](src/parallel_utils) | parallel utilities for accelerating computations of graph algorithms
[src/util](src/util) | utilities for graph input and pretty printing

## Scripts
content | detail
--- | ---
[scripts/analyze_algo_quality.py](scripts/analyze_algo_quality.py) | quality analyzer
[scripts/metrics/link_belong_modularity.py](scripts/metrics/link_belong_modularity.py) | link belonging modularity

## Demo Files
content | detail
--- | ---
[demo_files/demo_graph.csv](demo_files/demo_graph.csv) | toy graph in edge-list format
[demo_files/demo_result.txt](demo_files/demo_result.txt) | computation result of connected-iterative-scan algorithm

## Dataset
- [small_datasets/collaboration_edges_input.csv](small_datasets/collaboration_edges_input.csv), from uci repo
- [small_datasets/karate_edges_input.csv](small_datasets/karate_edges_input.csv), from uci repo
