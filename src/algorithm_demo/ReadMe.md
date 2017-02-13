#Algorithm Executable
##Algorithm Implementation Codes
- see [../algorithm](../algorithm)

##Main Entry Source Codes
- [demo_cis.cpp](demo_cis.cpp)
- [demo_demon.cpp](demo_demon.cpp)
- [demo_hk_growth.cpp](demo_hk_growth.cpp)

##Dependency
- boost regex module, requiring to link the executable to related run-time library

##Arguments Sample

- one argument specify the path of demo_graph, in edge list with weight representation

```zsh
/home/cheyulin/GitRepos/LocalityBasedGraphAlgo/demo_files/demo_graph.csv 
```

##Input File Sample

```zsh
#src_name dst_name edge_weight
0 1 3
1 2 1
2 3 4
3 4 4
0 5 3
2 7 4
4 7 4
4 8 1
```

##Output Sample

- cis algorithm

```zsh
#src_name dst_name edge_weight
0 1 3
1 2 1
2 3 4
3 4 4
0 5 3
2 7 4
4 7 4
4 8 1
src:0,dst:1,weight:3
src:1,dst:2,weight:1
src:2,dst:3,weight:4
src:3,dst:4,weight:4
src:0,dst:5,weight:3
src:2,dst:7,weight:4
src:4,dst:7,weight:4
src:4,dst:8,weight:1
idx result:[[0, 1, 5], [1, 2, 3, 4, 6, 7]]
name result:[[0, 1, 5], [1, 2, 3, 4, 7, 8]]
comm size:2
```
