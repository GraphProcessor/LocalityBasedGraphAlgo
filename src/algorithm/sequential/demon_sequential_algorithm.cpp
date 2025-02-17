//
// Created by cheyulin on 12/22/16.
//

#include "demon_sequential_algorithm.h"

namespace yche {
    unique_ptr<Demon::SubGraph> Demon::ExtractEgoMinusEgo(const Demon::Vertex &ego_vertex) const {
        auto ego_net_ptr = make_unique<SubGraph>();
        auto sub_vertex_id_map = get(vertex_id, *ego_net_ptr);
        auto sub_vertex_index_map = get(vertex_index, *ego_net_ptr);
        auto sub_vertex_weight_map = get(vertex_weight, *ego_net_ptr);
        auto sub_vertex_label_map = get(vertex_label, *ego_net_ptr);
        auto vertex_index_map = get(vertex_index, *graph_ptr_);
        auto vertex_weight_map = get(vertex_weight, *graph_ptr_);

        auto sub_vertices = vector<SubGraphVertex>();
        auto sub_graph_index_dict = map<int, int>();
        for (auto vp = adjacent_vertices(ego_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto sub_vertex = add_vertex(*ego_net_ptr);
            sub_vertices.emplace_back(sub_vertex);

            auto graph_vertex_index = static_cast<int>(vertex_index_map[*vp.first]);
            sub_vertex_id_map[sub_vertex] = graph_vertex_index;
            sub_vertex_weight_map[sub_vertex] = vertex_weight_map[*vp.first];
            sub_vertex_label_map[sub_vertex][0] = sub_vertex_id_map[sub_vertex];
            sub_vertex_label_map[sub_vertex][1] = 0;

            sub_graph_index_dict.emplace(graph_vertex_index, sub_vertex_index_map[sub_vertex]);
        }

        for (auto vp = adjacent_vertices(ego_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto source_sub_vertex_index = sub_graph_index_dict[vertex_index_map[*vp.first]];
            auto source_sub_vertex = sub_vertices[source_sub_vertex_index];
            for (auto vp_inner = adjacent_vertices(ego_vertex, *graph_ptr_);
                 vp_inner.first != vp_inner.second; ++vp_inner.first) {
                bool is_edge_exists = edge(*vp.first, *vp_inner.first, *graph_ptr_).second;
                if (is_edge_exists) {
                    auto end_sub_vertex_index = sub_graph_index_dict[vertex_index_map[*vp_inner.first]];
                    auto end_vertex = sub_vertices[end_sub_vertex_index];
                    add_edge(source_sub_vertex, end_vertex, *ego_net_ptr);
                }
            }
        }
        return ego_net_ptr;
    }

    void Demon::PropagateLabelSingle(const unique_ptr<SubGraph> &sub_graph_ptr, const SubGraphVertex &sub_graph_vertex,
                                     mt19937 &rand_generator, int last_label_idx, int curr_label_idx,
                                     property_map<SubGraph, vertex_weight_t>::type &sub_vertex_weight_map,
                                     property_map<SubGraph, vertex_label_t>::type &sub_vertex_label_map) const {
        auto label_weight_map = map<int, double>();
        //Label Propagation
        for (auto vp_inner = adjacent_vertices(sub_graph_vertex, *sub_graph_ptr);
             vp_inner.first != vp_inner.second; ++vp_inner.first) {
            auto neighbor_vertex = *vp_inner.first;
            auto neighbor_vertex_label = sub_vertex_label_map[neighbor_vertex][last_label_idx];
            auto neighbor_vertex_weight = sub_vertex_weight_map[neighbor_vertex];

            auto my_iterator = label_weight_map.find(neighbor_vertex_label);
            if (my_iterator == label_weight_map.end()) {
                label_weight_map.emplace(neighbor_vertex_label, neighbor_vertex_weight);
            } else {
                label_weight_map[neighbor_vertex_label] += neighbor_vertex_weight;
            }
        }

        //Find Maximum Vote
        auto candidate_label_vec = vector<int>();
        auto max_val = 0.0;
        if (label_weight_map.size() == 0) {
            sub_vertex_label_map[sub_graph_vertex][curr_label_idx] = sub_vertex_label_map[sub_graph_vertex][last_label_idx];
        } else {
            for (auto label_to_weight_pair:label_weight_map) {
                auto label_weight = label_to_weight_pair.second;
                if (label_weight > max_val) {
                    candidate_label_vec.clear();
                    max_val = label_weight;
                }
                if (label_weight >= max_val) {
                    candidate_label_vec.emplace_back(label_to_weight_pair.first);
                }
            }

            auto choice_index = 0;
            if (candidate_label_vec.size() >= 1) {
                auto distribution = uniform_int_distribution<>(0, static_cast<int>((candidate_label_vec.size() - 1)));
                choice_index = distribution(rand_generator);
            }
            sub_vertex_label_map[sub_graph_vertex][curr_label_idx] = candidate_label_vec[choice_index];
        }
    }

    Demon::CommunityVec Demon::GetCommunityVec(const unique_ptr<SubGraph> &sub_graph_ptr, const Vertex &ego_vertex,
                                               int curr_label_idx) const {
        auto sub_vertex_label_map = get(vertex_label, *sub_graph_ptr);
        auto sub_vertex_id_map = get(vertex_id, *sub_graph_ptr);
        auto vertex_index_map = get(vertex_index, *graph_ptr_);
        auto community_dict = map<int, Community>();

        for (auto vp = vertices(*sub_graph_ptr); vp.first != vp.second; ++vp.first) {
            auto sub_vertex = *vp.first;
            auto v_label = sub_vertex_label_map[sub_vertex][curr_label_idx];
            if (community_dict.find(v_label) == community_dict.end()) {
                community_dict.emplace(v_label, Community());
            }
            community_dict[v_label].emplace_back(sub_vertex_id_map[sub_vertex]);
        }

        auto community_vec = CommunityVec();
        for (auto &idx_community:community_dict) {
            //Add Ego Vertex
            auto &community = idx_community.second;
            community.emplace_back(vertex_index_map[ego_vertex]);
            sort(community.begin(), community.end());
            community_vec.emplace_back(std::move(community));
        }

        //Outlier
        if (community_dict.size() == 0) {
            community_vec.emplace_back(Community(1, static_cast<int>(vertex_index_map[ego_vertex])));
        }
        return community_vec;
    }

    Demon::CommunityVec Demon::PropagateLabel(const unique_ptr<SubGraph> &sub_graph_ptr,
                                              const Vertex &ego_vertex) const {
        auto sub_vertex_weight_map = get(vertex_weight, *sub_graph_ptr);
        auto sub_vertex_label_map = get(vertex_label, *sub_graph_ptr);
        static thread_local random_device rand_d;
        static thread_local std::mt19937 rand_generator(rand_d());

        auto last_label_idx = 0;
        auto curr_label_idx = 1;

        for (auto iter_num = 0; iter_num < max_iter_; iter_num++) {
            auto sub_vertices = vector<SubGraphVertex>();
            for (auto vp = vertices(*sub_graph_ptr); vp.first != vp.second; ++vp.first) {
                sub_vertices.emplace_back(*vp.first);
            }
            shuffle(sub_vertices.begin(), sub_vertices.end(), rand_generator);
            for (auto &sub_vertex:sub_vertices) {
                PropagateLabelSingle(sub_graph_ptr, sub_vertex, rand_generator, last_label_idx,
                                     curr_label_idx, sub_vertex_weight_map, sub_vertex_label_map);
            }
            swap(last_label_idx, curr_label_idx);
        }

        auto community_vec = GetCommunityVec(sub_graph_ptr, ego_vertex, curr_label_idx);
        return community_vec;
    }

    double Demon::GetIntersectRatio(const Community &left_community, const Community &right_community) {
        auto intersect_set = vector<int>(left_community.size() + right_community.size());
        auto iter_end = set_intersection(left_community.begin(), left_community.end(), right_community.begin(),
                                         right_community.end(), intersect_set.begin());
        auto inter_set_size = iter_end - intersect_set.begin();
        auto rate = static_cast<double>(inter_set_size) / min(left_community.size(), right_community.size());
        return rate;
    }

    Demon::Community Demon::GetUnion(const Community &left_community, const Community &right_community) {
        auto union_set = vector<int>(left_community.size() + right_community.size());
        auto iter_end = set_union(left_community.begin(), left_community.end(), right_community.begin(),
                                  right_community.end(), union_set.begin());
        union_set.resize(iter_end - union_set.begin());
        return union_set;
    }

    void Demon::MergeToGlobal(CommunityVec &result) {
        if (overlap_community_vec_.size() == 0) {
            for (auto &community:result) {
                if (community.size() > min_comm_size_) { overlap_community_vec_.emplace_back(std::move(community)); }
            }
        } else {
            for (auto &new_community:result) {
                auto first_access_flag = false;
                for (auto &old_community:overlap_community_vec_) {
                    if (GetIntersectRatio(new_community, old_community) > epsilon_) {
                        old_community = GetUnion(new_community, old_community);
                        break;
                    } else if (new_community.size() > min_comm_size_ && !first_access_flag) {
                        first_access_flag = true;
                    }
                }
                if (first_access_flag) { overlap_community_vec_.emplace_back(std::move(new_community)); }
            }
        }
    }

    Demon::CommunityVec Demon::ExecuteDemon() {
        for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto ego_vertex = *vp.first;
            auto sub_graph_ptr = ExtractEgoMinusEgo(ego_vertex);
            auto community_vec = PropagateLabel(sub_graph_ptr, ego_vertex);
            MergeToGlobal(community_vec);
        }
        return overlap_community_vec_;
    }
}