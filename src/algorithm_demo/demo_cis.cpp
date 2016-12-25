//
// Created by cheyulin on 10/5/16.
//

#include "algorithm/sequential/cis_sequential_algorithm.h"

#include "util/pretty_print.h"
#include "util/graph_io_helper.h"
#include "util/basic_io_helper.h"

using yche::Cis;
using yche::Edge;
using Vertex=Cis::Vertex;
using Graph=Cis::Graph;
using namespace std;

unique_ptr<Graph> ConstructGraph(unordered_map<int, Vertex> &vertex_dict,
                                 unordered_map<int, int> &name_dict, vector<Edge> &edges_vec) {
    auto graph_ptr = make_unique<Graph>();
    auto edge_weight_map = boost::get(boost::edge_weight, *graph_ptr);
    for (auto &my_edge : edges_vec) {
        if (vertex_dict.find(my_edge.src_index_) == vertex_dict.end()) {
            Vertex vertex = boost::add_vertex(*graph_ptr);
            vertex_dict.emplace(my_edge.src_index_, vertex);
        }
        if (vertex_dict.find(my_edge.dst_index_) == vertex_dict.end()) {
            Vertex vertex = boost::add_vertex(*graph_ptr);
            vertex_dict.emplace(my_edge.dst_index_, vertex);
        }
        auto edge = Cis::Edge();
        auto flag = false;
        tie(edge, flag) = add_edge(vertex_dict[my_edge.src_index_], vertex_dict[my_edge.dst_index_], *graph_ptr);
        if (flag) {
            edge_weight_map[edge] = my_edge.edge_weight_;
            cout << my_edge;
        }
    }

    auto vertex_index_map = boost::get(boost::vertex_index, *graph_ptr);
    for (auto iter = vertex_dict.begin(); iter != vertex_dict.end(); ++iter) {
        name_dict.emplace(vertex_index_map[iter->second], iter->first);
    }
    return graph_ptr;
}

int main(int argc, char *argv[]) {
    auto edges_vec = yche::ReadWeightedEdgeList(argv[1]);

    auto vertex_dict = unordered_map<int, Vertex>();
    auto name_dict = unordered_map<int, int>();
    auto graph_ptr = ConstructGraph(vertex_dict, name_dict, edges_vec);

    auto cis = Cis(graph_ptr, 0);
    cis.ExecuteCis();

    auto &arr_2d = cis.overlap_community_vec_;
    auto name_arr_2d = yche::Map2DArrWithDict(arr_2d, name_dict);
    cout << "idx result:" << arr_2d << endl;
    cout << "name result:" << name_arr_2d << endl;
    return 0;
}