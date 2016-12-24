//
// Created by cheyulin on 12/22/16.
//

#include "algorithm/sequential/demon_sequential_algorithm.h"
#include "util/graph_io_helper.h"
#include "util/pretty_print.h"

using namespace std;
using yche::Demon;
using Graph=Demon::Graph;
using Vertex=Demon::Vertex;

void ConstructGraphWithEdgeVecForDemon(unique_ptr<Graph> &graph_ptr, map<int, Vertex> &name_vertex_map,
                                       map<int, int> &index_name_map, vector<pair<int, int>> &edges_vec) {
    auto vertex_weight_map = boost::get(boost::vertex_weight, *graph_ptr);
    for (auto iter = edges_vec.begin(); iter != edges_vec.end(); ++iter) {
        if (name_vertex_map.find(iter->first) == name_vertex_map.end()) {
            Vertex vertex = add_vertex(*graph_ptr);
            vertex_weight_map[vertex] = 1;
            name_vertex_map.insert(make_pair(iter->first, vertex));
        }
        if (name_vertex_map.find(iter->second) == name_vertex_map.end()) {
            Vertex vertex = add_vertex(*graph_ptr);
            vertex_weight_map[vertex] = 1;
            name_vertex_map.insert(make_pair(iter->second, vertex));
        }
        add_edge(name_vertex_map[iter->first], name_vertex_map[iter->second], *graph_ptr);
    }

    auto vertex_index_map = boost::get(boost::vertex_index, *graph_ptr);
    for (auto iter = name_vertex_map.begin(); iter != name_vertex_map.end(); ++iter) {
        index_name_map.insert(make_pair(vertex_index_map[iter->second], iter->first));
    }
}

int main(int argc, char *argv[]) {
    auto edges_vec = yche::ReadEdgeList(argv[1]);

    auto graph_ptr = make_unique<Graph>();
    map<int, Vertex> name_vertex_map;
    map<int, int> index_name_map;
    ConstructGraphWithEdgeVecForDemon(graph_ptr, name_vertex_map, index_name_map, edges_vec);

    auto epsilon = 0.25;
    auto min_community_size = 3;
    auto max_iteration = 100;
    auto demon_algo = Demon(epsilon, min_community_size, graph_ptr, max_iteration);

    demon_algo.ExecuteDaemon();
    cout << demon_algo.overlap_community_vec_ << endl;
}