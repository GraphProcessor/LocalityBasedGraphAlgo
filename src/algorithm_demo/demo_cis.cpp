//
// Created by cheyulin on 10/5/16.
//

#include "algorithm/cis_sequential_algorithm.h"
#include "input_output_handler.h"
#include "util/pretty_print.h"

using namespace yche;

void ConstructGraphWithEdgeVecForCIS(unique_ptr<Cis::Graph> &graph_ptr, map<int, Cis::Vertex> &name_vertex_map,
                                     map<int, int> &index_name_map, vector<Edge> &edges_vec) {
    using namespace boost;
    property_map<Cis::Graph, edge_weight_t>::type edge_weight_map = get(edge_weight, *graph_ptr);
    cout << endl;
    for (auto iter = edges_vec.begin(); iter != edges_vec.end(); ++iter) {
        if (name_vertex_map.find(iter->src_index_) == name_vertex_map.end()) {
            Cis::Vertex vertex = add_vertex(*graph_ptr);
            name_vertex_map.insert(make_pair(iter->src_index_, vertex));
        }
        if (name_vertex_map.find(iter->dst_index_) == name_vertex_map.end()) {
            Cis::Vertex vertex = add_vertex(*graph_ptr);
            name_vertex_map.insert(make_pair(iter->dst_index_, vertex));
        }
        Cis::Edge edge;
        bool flag;
        tie(edge, flag) = add_edge(name_vertex_map[iter->src_index_], name_vertex_map[iter->dst_index_], *graph_ptr);
        edge_weight_map[edge] = iter->edge_weight_;
        cout << "src:" << iter->src_index_ << ",dst:" << iter->dst_index_ << ",weight:" << iter->edge_weight_ << endl;
    }

    property_map<Cis::Graph, vertex_index_t>::type vertex_index_map = get(vertex_index, *graph_ptr);
    for (auto iter = name_vertex_map.begin(); iter != name_vertex_map.end(); ++iter) {
        index_name_map.insert(make_pair(vertex_index_map[iter->second], iter->first));
    }
}

int main(int argc, char *argv[]) {
    auto file_name = argv[1];
    vector<Edge> edges_vec;
    ReadEdgeListWithWeightInToEdgeVector(file_name, edges_vec);

    auto graph_ptr = make_unique<Cis::Graph>();
    map<int, Cis::Vertex> name_vertex_map;
    map<int, int> index_name_map;
    ConstructGraphWithEdgeVecForCIS(graph_ptr, name_vertex_map, index_name_map, edges_vec);

    auto cis = Cis(graph_ptr, 0);
    cis.ExecuteCis();
    cout << cis.overlap_community_vec_ << endl;
    return 0;
}