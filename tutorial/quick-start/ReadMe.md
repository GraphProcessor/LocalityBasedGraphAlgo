#Quick Start
##Recommended IDEs

Jetbrains IDEs are easy to learn and flexible to use. You can use 
[jetbrains-toolbox](https://www.jetbrains.com/toolbox/app/?fromMenu) to help you manage the update of following 
useful ides.

tool | language
--- | ---
[clion](https://www.jetbrains.com/clion/?fromMenu)  | for c++
[pycharm](https://www.jetbrains.com/pycharm/?fromMenu)  | for python
[webstorm](https://www.jetbrains.com/webstorm/?fromMenu) | for javascript


##Build Tool Chain

tool | functionality | links
--- | --- | ---
gcc4.9+ | compile source codes with std=c++14 | 
make | build codes with makefile | 
cmake | generate makefile | [cmake-tutorial](https://cmake.org/cmake-tutorial/)

##Dependency
- boost graph library, only requiring including headers, which use template-meta-programming. Adjacency list 
and compressed sparse row graph concept are used in this project.

configure as following in [CMakeLists.txt](../../src/CMakeLists.txt) to give gcc 
additional-search-path for header-files.

```
include_directories(${Boost_INCLUDE_DIRS})
```

- boost regex, requiring link of the related dynamic link library

used in [../src/util/graph_io_helper.h](../src/util/graph_io_helper.h) for graph input handling.

##Build Steps(On Linux), with boost-dev and build-tool-chain installed

in root directory of this project

- first build the project to have executables

```zsh
mkdir -p build
cd build
cmake ../src/
make -j
```

- second run some executables on datasets

e.g, cis algorithm 

enter build directory, then do as follows, where the first argument is `../small_datasets/collaboration_edges_input.csv`

```zsh
algorithm_demo/demo_cis ../small_datasets/collaboration_edges_input.csv
```