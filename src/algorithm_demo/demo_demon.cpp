//
// Created by cheyulin on 12/22/16.
//

#include "algorithm/sequential/demon_sequential_algorithm.h"
#include "util/graph_io_helper.h"
#include "util/basic_io_helper.h"
#include "util/pretty_print.h"

using namespace std;
using yche::Demon;
using Graph=Demon::Graph;
using Vertex=Demon::Vertex;

unique_ptr<Graph> ConstructGraph(map<int, Vertex> &name_vertex_map, map<int, int> &index_name_map,
                                 vector<pair<int, int>> &edges_vec) {
    auto graph_ptr = make_unique<Graph>();
    auto vertex_weight_map = boost::get(boost::vertex_weight, *graph_ptr);
    for (auto &edge:edges_vec) {
        if (name_vertex_map.find(edge.first) == name_vertex_map.end()) {
            Vertex vertex = add_vertex(*graph_ptr);
            vertex_weight_map[vertex] = 1;
            name_vertex_map.insert(make_pair(edge.first, vertex));
        }
        if (name_vertex_map.find(edge.second) == name_vertex_map.end()) {
            Vertex vertex = add_vertex(*graph_ptr);
            vertex_weight_map[vertex] = 1;
            name_vertex_map.insert(make_pair(edge.second, vertex));
        }
        auto edge_flag_pair = add_edge(name_vertex_map[edge.first], name_vertex_map[edge.second], *graph_ptr);
        if (edge_flag_pair.second) {
            cout << "edge:" << edge_flag_pair.first << endl;
        }
    }

    auto vertex_index_map = boost::get(boost::vertex_index, *graph_ptr);
    for (auto &name_vertex:name_vertex_map) {
        index_name_map.emplace(vertex_index_map[name_vertex.second], name_vertex.first);
    }
    return graph_ptr;
}

int main(int argc, char *argv[]) {
    auto edges_vec = yche::ReadEdgeList(argv[1]);

    auto vertex_dict = map<int, Vertex>();
    auto name_dict = map<int, int>();

    auto epsilon = 0.25;
    auto min_community_size = 1;
    auto max_iteration = 100;
    auto demon_algo = Demon(epsilon, min_community_size, ConstructGraph(vertex_dict, name_dict, edges_vec),
                            max_iteration);

    auto arr_2d = std::move(demon_algo.ExecuteDemon());
    auto name_arr_2d = yche::Map2DArrWithDict(arr_2d, name_dict);

    cout << "idx result:" << arr_2d << endl;
    cout << "name result:" << name_arr_2d << endl;
    cout << "comm size:" << name_arr_2d.size() << endl;
}