//
// Created by cheyulin on 12/22/16.
//

#ifndef CODES_YCHE_DEMON_SEQUENTIAL_ALGORITHM_H
#define CODES_YCHE_DEMON_SEQUENTIAL_ALGORITHM_H

#include <memory>
#include <vector>
#include <random>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>

#include "util/pretty_print.h"

namespace boost {
    enum vertex_weight_t {
        vertex_weight
    };
    enum vertex_label_t {
        vertex_label
    };
    enum vertex_id_t {
        vertex_id
    };

    BOOST_INSTALL_PROPERTY(vertex, weight);
    BOOST_INSTALL_PROPERTY(vertex, label);
    BOOST_INSTALL_PROPERTY(vertex, id);
}

namespace yche {
    using namespace boost;
    using namespace std;

    class Demon {
    public:
        using VertexProperties = property<vertex_weight_t, double, property<vertex_index_t, int>>;
        using Graph = adjacency_list<hash_setS, vecS, undirectedS, VertexProperties>;
        using Vertex = graph_traits<Graph>::vertex_descriptor;

        using SubGraphVertexProperties = property<vertex_weight_t, double,
                property<vertex_id_t, int, property<vertex_label_t, array<int, 2>>>>;
        using SubGraph = adjacency_list<hash_setS, vecS, undirectedS, SubGraphVertexProperties>;
        using SubGraphVertex = graph_traits<SubGraph>::vertex_descriptor;

        using Community = vector<int>;
        using CommunityVec = vector<Community>;

        CommunityVec overlap_community_vec_;

        Demon(double epsilon, int min_comm_size, unique_ptr<Graph> &graph_ptr, int max_iter) :
                epsilon_(epsilon), min_comm_size_(min_comm_size), max_iter_(max_iter),
                graph_ptr_(std::move(graph_ptr)) {}

        void ExecuteDemon();

    private:
        unique_ptr<Graph> graph_ptr_;
        double epsilon_;
        int min_comm_size_;
        int max_iter_;

        unique_ptr<SubGraph> ExtractEgoMinusEgo(Vertex &ego_vertex) const;

        void PropagateLabelSingle(unique_ptr<SubGraph> &sub_graph_ptr, SubGraphVertex &sub_graph_vertex,
                                  mt19937 &rand_generator, int last_label_idx, int curr_label_idx,
                                  property_map<SubGraph, vertex_weight_t>::type &sub_vertex_weight_map,
                                  property_map<SubGraph, vertex_label_t>::type &sub_vertex_label_map);

        CommunityVec GetCommunityVec(unique_ptr<SubGraph> &sub_graph_ptr, Vertex &ego_vertex, int curr_label_idx);

        CommunityVec PropagateLabel(unique_ptr<SubGraph> &sub_graph_ptr, Vertex &ego_vertex);

        double GetIntersectRatio(Community &left_community, Community &right_community) const;

        Community GetUnion(Community &left_community, Community &right_community) const;

        void MergeToGlobal(CommunityVec &result);
    };
}

#endif //CODES_YCHE_DEMON_SEQUENTIAL_ALGORITHM_H
