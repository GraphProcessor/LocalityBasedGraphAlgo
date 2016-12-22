//
// Created by cheyulin on 12/22/16.
//

#include "demon_sequential_algorithm.h"

namespace yche {
    unique_ptr<Demon::SubGraph> Demon::ExtractEgoMinusEgo(Demon::Vertex &ego_vertex) const {
        auto ego_net_ptr = make_unique<SubGraph>();
        property_map<SubGraph, vertex_id_t>::type sub_vertex_id_map = get(vertex_id, *ego_net_ptr);
        property_map<SubGraph, vertex_index_t>::type sub_vertex_index_map = get(vertex_index, *ego_net_ptr);
        property_map<SubGraph, vertex_weight_t>::type sub_vertex_weight_map = get(vertex_weight, *ego_net_ptr);
        property_map<SubGraph, vertex_label_t>::type sub_vertex_label_map = get(vertex_label, *ego_net_ptr);
        property_map<Graph, vertex_index_t>::type vertex_index_map = get(vertex_index, *graph_ptr_);
        property_map<Graph, vertex_weight_t>::type vertex_weight_map = get(vertex_weight, *graph_ptr_);

        auto sub_vertices = vector<SubGraphVertex>();
        auto sub_graph_index_dict = map<int, int>();
        for (auto vp = adjacent_vertices(ego_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            graph_traits<SubGraph>::vertex_descriptor sub_vertex = add_vertex(*ego_net_ptr);
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
                    auto end_vertex_ptr = sub_vertices[end_sub_vertex_index];
                    add_edge(source_sub_vertex, end_vertex_ptr, *ego_net_ptr);
                }
            }
        }
        return ego_net_ptr;
    }

    void Demon::PropagateLabelSingle(unique_ptr<SubGraph> &sub_graph_ptr, SubGraphVertex &sub_graph_Vertex,
                                     std::mt19937 &rand_generator, int last_index_indicator, int curr_index_indicator,
                                     property_map<SubGraph, vertex_weight_t>::type &sub_vertex_weight_map,
                                     property_map<SubGraph, vertex_label_t>::type &sub_vertex_label_map) {
        auto label_weight_map = map<int, double>();
        //Label Propagation
        for (auto vp_inner = adjacent_vertices(sub_graph_Vertex, *sub_graph_ptr);
             vp_inner.first != vp_inner.second; ++vp_inner.first) {
            auto neighbor_vertex = *vp_inner.first;
            auto neighbor_vertex_label = sub_vertex_label_map[neighbor_vertex][last_index_indicator];
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
        auto current_vertex = sub_graph_Vertex;
        if (label_weight_map.size() == 0) {
            sub_vertex_label_map[current_vertex][curr_index_indicator] = sub_vertex_label_map[current_vertex][last_index_indicator];
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
            sub_vertex_label_map[current_vertex][curr_index_indicator] = candidate_label_vec[choice_index];
        }
    }

    Demon::CommunityVec Demon::GetCommunityVec(unique_ptr<SubGraph> &sub_graph_ptr,
                                               Vertex &ego_vertex, int curr_index_indicator) {
        property_map<SubGraph, vertex_label_t>::type sub_vertex_label_map = get(vertex_label, *sub_graph_ptr);
        property_map<SubGraph, vertex_id_t>::type sub_vertex_id_map = get(vertex_id, *sub_graph_ptr);
        property_map<Graph, vertex_index_t>::type vertex_index_map = get(vertex_index, *graph_ptr_);

        map<int, Community> label_indices_map;
        for (auto vp = vertices(*sub_graph_ptr); vp.first != vp.second; ++vp.first) {
            auto sub_vertex = *vp.first;
            auto v_label = sub_vertex_label_map[sub_vertex][curr_index_indicator];
            if (label_indices_map.find(v_label) == label_indices_map.end()) {
                label_indices_map.emplace(v_label, Community());
            }
            label_indices_map[v_label].emplace_back(sub_vertex_id_map[sub_vertex]);
        }

        CommunityVec community_vec = CommunityVec();
        for (auto iter = label_indices_map.begin(); iter != label_indices_map.end(); ++iter) {
            //Add Ego Vertex
            //Make The Community Vector Sorted
            iter->second.emplace_back(vertex_index_map[ego_vertex]);
            sort(iter->second.begin(), iter->second.end());
            community_vec.emplace_back(std::move(iter->second));
        }
        if (label_indices_map.size() == 0) {
            //Outlier
            Community comm_ptr = Community();
            comm_ptr.emplace_back(vertex_index_map[ego_vertex]);
            community_vec.emplace_back(std::move(comm_ptr));
        }
        return community_vec;
    }

    Demon::CommunityVec Demon::PropagateLabel(unique_ptr<SubGraph> &sub_graph_ptr, Vertex &ego_vertex) {
        property_map<SubGraph, vertex_weight_t>::type sub_vertex_weight_map = get(vertex_weight, *sub_graph_ptr);
        property_map<SubGraph, vertex_label_t>::type sub_vertex_label_map = get(vertex_label, *sub_graph_ptr);
        auto iteration_num = 0;
        for (; iteration_num < max_iter_; iteration_num++) {
            auto curr_index_indicator = (iteration_num + 1) % 2;
            auto last_index_indicator = iteration_num % 2;
            auto all_sub_vertices = vector<SubGraphVertex>();
            for (auto vp = vertices(*sub_graph_ptr); vp.first != vp.second; ++vp.first) {
                all_sub_vertices.emplace_back(*vp.first);
            }
            static thread_local random_device rand_d;
            static thread_local std::mt19937 rand_generator(rand_d());
            shuffle(all_sub_vertices.begin(), all_sub_vertices.end(), rand_generator);

            for (auto vertex_iter = all_sub_vertices.begin(); vertex_iter != all_sub_vertices.end(); ++vertex_iter) {
                PropagateLabelSingle(sub_graph_ptr, *vertex_iter, rand_generator, last_index_indicator,
                                     curr_index_indicator, sub_vertex_weight_map, sub_vertex_label_map);
            }
        }

        auto curr_index_indicator = (iteration_num + 1) % 2;
        auto community_vec = GetCommunityVec(sub_graph_ptr, ego_vertex, curr_index_indicator);
        return community_vec;
    }

    double Demon::GetIntersectRatio(Community &left_community, Community &right_community) {
        vector<int> intersect_set(left_community.size() + right_community.size());
        auto iter_end = set_intersection(left_community.begin(), left_community.end(), right_community.begin(),
                                         right_community.end(), intersect_set.begin());
        auto inter_set_size = iter_end - intersect_set.begin();
        double rate = static_cast<double>(inter_set_size) / min(left_community.size(), right_community.size());
        return rate;
    }

    Community Demon::GetUnion(Community &left_community, Community &right_community) {
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
                bool first_access_flag = false;
                for (auto &old_community:overlap_community_vec_) {
                    auto cover_rate_result = GetIntersectRatio(new_community, old_community);
                    if (cover_rate_result > epsilon_) {
                        old_community = GetUnion(new_community, old_community);
                        break;
                    } else if (new_community.size() > min_comm_size_ && !first_access_flag) {
                        first_access_flag = true;
                    }
                }
                if (first_access_flag) {
                    overlap_community_vec_.emplace_back(std::move(new_community));
                }
            }
        }
    }

    void Demon::ExecuteDaemon() {
        for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto ego_vertex = *vp.first;
            auto sub_graph_ptr = ExtractEgoMinusEgo(ego_vertex);
            auto community_vec = PropagateLabel(std::move(sub_graph_ptr), ego_vertex);
            MergeToGlobal(community_vec);
        }
    }
}