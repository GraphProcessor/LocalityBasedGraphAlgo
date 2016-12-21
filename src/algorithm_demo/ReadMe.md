##Arguments Sample

- one argument specify the path of demo_graph, in edge list with weight representation

```zsh
/home/cheyulin/GitRepos/LocalityBasedGraphAlgo/demo_files/demo_graph.csv 
```

##Debug Info

- error new one

```zsh
0 1 3
src:0,dst:1,edge_weight:3
1 2 1
src:1,dst:2,edge_weight:1
2 3 4
src:2,dst:3,edge_weight:4
3 4 4
src:3,dst:4,edge_weight:4
0 5 3
src:0,dst:5,edge_weight:3
2 7 4
src:2,dst:7,edge_weight:4
4 7 4
src:4,dst:7,edge_weight:4
4 8 1
src:4,dst:8,edge_weight:1
 src:0,dst:1,weight:3
 src:1,dst:2,weight:1
 src:2,dst:3,weight:4
 src:3,dst:4,weight:4
 src:0,dst:5,weight:3
 src:2,dst:7,weight:4
 src:4,dst:7,weight:4
 src:4,dst:8,weight:1
[[0, 1, 2, 3, 4, 5, 6, 7]]
```

- correct previous one 

```zsh
0 1 3
src:0,dst:1,edge_weight:3
1 2 1
src:1,dst:2,edge_weight:1
2 3 4
src:2,dst:3,edge_weight:4
3 4 4
src:3,dst:4,edge_weight:4
0 5 3
src:0,dst:5,edge_weight:3
2 7 4
src:2,dst:7,edge_weight:4
4 7 4
src:4,dst:7,edge_weight:4
4 8 1
src:4,dst:8,edge_weight:1
 src:0,dst:1,weight:3
 src:1,dst:2,weight:1
 src:2,dst:3,weight:4
 src:3,dst:4,weight:4
 src:0,dst:5,weight:3
 src:2,dst:7,weight:4
 src:4,dst:7,weight:4
 src:4,dst:8,weight:1
0(0),1(1),5(5),
1(1),2(2),3(3),4(4),7(6),8(7),
```