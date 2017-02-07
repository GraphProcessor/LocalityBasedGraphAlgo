//
// Created by cheyulin on 1/5/17.
//

#include "algorithm/sequential/hk_grow_sequential_algorithm.h"

#include "util/graph_io_helper.h"
#include "util/basic_io_helper.h"
#include "util/pretty_print.h"

using yche::HKGrow;
using Graph=HKGrow::Graph;
using Vertex=size_t;
using namespace std;

unique_ptr<Graph> ConstructGraph(map<int, Vertex> &name_vertex_map, map<Vertex, int> &index_name_map,
                                 vector<pair<int, int>> &edges_vec) {
    auto v_set = set<int>();
    for (auto &edge:edges_vec) {
        v_set.emplace(edge.first);
        v_set.emplace(edge.second);
    }

    auto iter_num = 0;
    for (auto name:v_set) {
        name_vertex_map.emplace(name, iter_num);
        index_name_map.emplace(iter_num, name);
        iter_num++;
    }

    auto edges = vector<pair<int, int>>();
    for (auto &edge:edges_vec) {
        edges.emplace_back(name_vertex_map[edge.first], name_vertex_map[edge.second]);
        edges.emplace_back(name_vertex_map[edge.second], name_vertex_map[edge.first]);
    }

    auto last = unique(begin(edges), end(edges));
    edges.erase(last, end(edges));
    cout << edges << endl;
    return make_unique<Graph>(boost::edges_are_unsorted, begin(edges), end(edges), v_set.size());
}

int main(int argc, char *argv[]) {
    auto edges_vec = yche::ReadEdgeList(argv[1]);
    auto vertex_dict = map<int, Vertex>();
    auto name_dict = map<Vertex, int>();

    auto hkgrow_algo = HKGrow(ConstructGraph(vertex_dict, name_dict, edges_vec), 15, 0.0001);
    auto arr_2d = std::move(hkgrow_algo.ExecuteHRGRow());
    cout << arr_2d << endl;
}