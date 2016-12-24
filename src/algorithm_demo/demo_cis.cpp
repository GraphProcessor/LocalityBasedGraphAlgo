//
// Created by cheyulin on 10/5/16.
//

#include "algorithm/sequential/cis_sequential_algorithm.h"
#include "input_output_handler.h"
#include "util/pretty_print.h"

using namespace yche;

using Vertex=Cis::Vertex;
using Graph=Cis::Graph;

void ConstructGraph(unique_ptr<Graph> &graph_ptr, std::unordered_map<int, Vertex> &vertex_dict,
                    std::unordered_map<int, int> &name_dict, vector<Edge> &edges_vec) {
    auto edge_weight_map = get(edge_weight, *graph_ptr);
    cout << '\n';
    for (auto &my_edge : edges_vec) {
        if (vertex_dict.find(my_edge.src_index_) == vertex_dict.end()) {
            Vertex vertex = add_vertex(*graph_ptr);
            vertex_dict.emplace(my_edge.src_index_, vertex);
        }
        if (vertex_dict.find(my_edge.dst_index_) == vertex_dict.end()) {
            Vertex vertex = add_vertex(*graph_ptr);
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

    auto vertex_index_map = get(vertex_index, *graph_ptr);
    for (auto iter = vertex_dict.begin(); iter != vertex_dict.end(); ++iter) {
        name_dict.emplace(vertex_index_map[iter->second], iter->first);
    }
}

int main(int argc, char *argv[]) {
    auto edges_vec = ReadWeightedEdgeList(argv[1]);

    auto graph_ptr = make_unique<Graph>();
    auto vertex_dict = std::unordered_map<int, Vertex>();
    auto name_dict = std::unordered_map<int, int>();
    ConstructGraph(graph_ptr, vertex_dict, name_dict, edges_vec);

    auto cis = Cis(graph_ptr, 0);
    cis.ExecuteCis();
    cout << cis.overlap_community_vec_ << endl;
    return 0;
}