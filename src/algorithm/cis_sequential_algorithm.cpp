//
// Created by cheyulin on 12/20/16.
//

#include "cis_sequential_algorithm.h"


namespace yche {
    Cis::Cis(unique_ptr<Cis::Graph> &graph_ptr, double lambda) : lambda_(lambda), graph_ptr_(std::move(graph_ptr)) {
        for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) { vertices_.emplace_back(*vp.first); }
    }

    double Cis::CalDensity(int size, double w_in, double w_out, double lambda) const {
        if (size < 1) {
            return numeric_limits<double>::min();
        } else {
            double partA = ((1 - lambda) * (w_in / (w_in + w_out)));
            double partB = (lambda * ((2 * w_in) / (size * (size - 1))));
            if (size == 1)
                partB = lambda;
            return partA + partB;
        }
    }

    double Cis::CalDensity(Community &community) const {
        return CalDensity(static_cast<int>(community.member_indices_.size()),
                          community.w_in_, community.w_out_, lambda_);
    }

    double Cis::CalDensity(Community &community, Member &member, MutationType mutation_type) const {
        if (mutation_type == MutationType::add_neighbor) {
            return CalDensity(static_cast<int>(community.member_indices_.size() + 1),
                              community.w_in_ + member.w_in_, community.w_out_ + member.w_out_, lambda_);

        } else {
            return CalDensity(static_cast<int>(community.member_indices_.size() - 1),
                              community.w_in_ - member.w_in_, community.w_out_ - member.w_out_, lambda_);
        }
    }

    void Cis::InitializeSeeds(const MemberSet &seed, Community &community,
                              MemberMap &member_dict, MemberMap &neighbor_dict,
                              MemberSet &to_computed_neighbors,
                              property_map<Graph, vertex_index_t>::type &vertex_index_map,
                              property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        for (auto &seed_vertex_index :seed) {
            auto member_info_ptr = Member(seed_vertex_index);
            community.member_indices_.emplace(seed_vertex_index);
            Vertex seed_vertex = vertices_[seed_vertex_index];

            for (auto vp = adjacent_vertices(seed_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                auto neighbor_vertex_index = static_cast<int>(vertex_index_map[*vp.first]);
                auto edge_weight = edge_weight_map[
                        edge(seed_vertex, vertices_[neighbor_vertex_index], *graph_ptr_).first];
                if (seed.find(neighbor_vertex_index) != seed.end()) {
                    member_info_ptr.w_in_ += edge_weight;
                    community.w_in_ += edge_weight;
                } else {
                    member_info_ptr.w_out_ += edge_weight;
                    community.w_out_ += edge_weight;
                    to_computed_neighbors.emplace(neighbor_vertex_index);
                }
            }
            member_dict.emplace(member_info_ptr.member_index_, member_info_ptr);
        }

        for (auto &neighbor_vertex_index :to_computed_neighbors) {
            auto neighbor_info_ptr = Member(neighbor_vertex_index);
            Vertex neighbor_vertex = vertices_[neighbor_vertex_index];
            for (auto vp = adjacent_vertices(neighbor_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                auto neighbor_neighbor_vertex_index = static_cast<int>(vertex_index_map[*vp.first]);
                auto edge_weight = edge_weight_map[edge(neighbor_vertex, vertices_[neighbor_neighbor_vertex_index],
                                                        *graph_ptr_).first];
                if (seed.find(neighbor_neighbor_vertex_index) != seed.end()) {
                    neighbor_info_ptr.w_in_ += edge_weight;
                } else {
                    neighbor_info_ptr.w_out_ += edge_weight;
                }
            }
            neighbor_dict.emplace(neighbor_info_ptr.member_index_, neighbor_info_ptr);
        }
    }

    void Cis::UpdateForAddNeighbor(const Cis::Vertex &mutate_vertex,
                                   Community &community,
                                   MemberMap &member_dict, MemberMap &neighbor_dict,
                                   property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                   property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        //Update Member and Neighbor List
        for (auto vp = adjacent_vertices(mutate_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto check_neighbor_vertex = *vp.first;
            auto check_neighbor_vertex_index = vertex_index_map[check_neighbor_vertex];
            auto check_neighbor_ptr = make_unique<Member>(check_neighbor_vertex_index);
            auto edge_weight = edge_weight_map[edge(mutate_vertex, check_neighbor_vertex,
                                                    *graph_ptr_).first];

            auto iter = member_dict.find(check_neighbor_ptr->member_index_);
            if (iter != member_dict.end() ||
                (iter = neighbor_dict.find(check_neighbor_ptr->member_index_)) != neighbor_dict.end()) {
                //Update Info In Members and Neighbors
                iter->second.w_in_ += edge_weight;
                iter->second.w_out_ -= edge_weight;
            } else {
                //Add New Neighbor
                auto member_info_ptr = Member(check_neighbor_vertex_index);
                for (auto vp_inner = adjacent_vertices(check_neighbor_vertex, *graph_ptr_);
                     vp_inner.first != vp_inner.second; ++vp_inner.first) {
                    auto neighbor_neighbor_vertex_index = static_cast<int>(vertex_index_map[*vp_inner.first]);
                    edge_weight = edge_weight_map[edge(check_neighbor_vertex,
                                                       vertices_[neighbor_neighbor_vertex_index], *graph_ptr_).first];
                    if (community.member_indices_.find(neighbor_neighbor_vertex_index) !=
                        community.member_indices_.end()) {
                        member_info_ptr.w_in_ += edge_weight;
                    } else {
                        member_info_ptr.w_out_ += edge_weight;
                    }
                }
                neighbor_dict.emplace(member_info_ptr.member_index_, member_info_ptr);
            }
        }
    }

    void Cis::UpdateForRemoveMember(const Cis::Vertex &mutate_vertex, Community &community,
                                    MemberMap &member_dict, MemberMap &neighbor_dict,
                                    property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                    property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        //Update Member and Neighbor List
        for (auto vp = adjacent_vertices(mutate_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto check_neighbor_vertex = *vp.first;
            auto check_neighbor_vertex_index = vertex_index_map[check_neighbor_vertex];
            auto check_neighbor_ptr = std::move(make_unique<Member>(check_neighbor_vertex_index));
            auto edge_weight = edge_weight_map[edge(mutate_vertex, check_neighbor_vertex, *graph_ptr_).first];

            auto iter = member_dict.find(check_neighbor_ptr->member_index_);
            if (iter != member_dict.end() ||
                (iter = neighbor_dict.find(check_neighbor_ptr->member_index_)) != neighbor_dict.end()) {
                //Update Info In Members and Neighbors
                iter->second.w_in_ -= edge_weight;
                iter->second.w_out_ += edge_weight;
            }
        }
    }

    void Cis::UpdateStates(const Cis::Vertex &mutate_vertex, Community &community,
                           MemberMap &member_dict, MemberMap &neighbor_dict, MutationType mutation_type,
                           property_map<Graph, vertex_index_t>::type &vertex_index_map,
                           property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        if (mutation_type == MutationType::add_neighbor) {
            UpdateForAddNeighbor(mutate_vertex, community, member_dict,
                                 neighbor_dict, vertex_index_map, edge_weight_map);
        } else {
            UpdateForRemoveMember(mutate_vertex, community, member_dict,
                                  neighbor_dict, vertex_index_map, edge_weight_map);
        }
    }


    Community Cis::SplitAndChoose(MemberSet &member_set) {
        vector<Community> community_vec;
        queue<int> frontier;
        std::unordered_set<int> mark_set;
        property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
        property_map<Graph, edge_weight_t>::type edge_weight_map = boost::get(edge_weight, *graph_ptr_);

        while (member_set.size() > 0) {
            int first_vertex = *member_set.begin();
            frontier.push(first_vertex);
            mark_set.emplace(first_vertex);

            auto community = Community();
            //One Connected Component
            while (frontier.size() > 0) {
                Edge edge;
                bool edge_exist_flag;
                auto expand_vertex_index = frontier.front();
                auto expand_vertex = vertices_[expand_vertex_index];
                //Add To Community Only At Start of Expansion
                community.member_indices_.emplace(expand_vertex_index);
                for (auto vp = adjacent_vertices(vertices_[expand_vertex_index], *graph_ptr_);
                     vp.first != vp.second; ++vp.first) {
                    auto neighbor_vertex = *vp.first;
                    auto adjacency_vertex_index = static_cast<int>(vertex_index_map[neighbor_vertex]);
                    tie(edge, edge_exist_flag) = boost::edge(expand_vertex, neighbor_vertex, *graph_ptr_);
                    //Do W_in and W_out Computation
                    auto iter = member_set.find(adjacency_vertex_index);
                    if (mark_set.find(adjacency_vertex_index) == mark_set.end() && iter != member_set.end()) {
                        community.w_in_ += edge_weight_map[edge];
                        //Erase and Enqueue
                        frontier.push(adjacency_vertex_index);
                        mark_set.emplace(adjacency_vertex_index);
                    } else {
                        community.w_out_ = edge_weight_map[edge];
                    }
                }
                member_set.erase(expand_vertex_index);
                frontier.pop();
            }
            community_vec.emplace_back(community);
        }

        sort(community_vec.begin(), community_vec.end(),
             [this](auto &left_comm, auto &right_comm) {
                 double left_density = this->CalDensity(left_comm);
                 double right_density = this->CalDensity(right_comm);
                 if (left_density != right_density) {
                     return left_density > right_density;
                 } else {
                     return (left_comm.member_indices_).size() > (right_comm.member_indices_).size();
                 }
             });
        return community_vec[0];
    }

    MemberVec Cis::ExpandSeed(MemberSet &seed) {
        //First: InitializeSeeds members and neighbors(or frontier)
        auto community = Community();
        MemberMap member_dict;
        MemberMap neighbor_dict;
        MemberSet to_computed_neighbors;
        property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
        property_map<Graph, edge_weight_t>::type edge_weight_map = boost::get(edge_weight, *graph_ptr_);

        InitializeSeeds(seed, community, member_dict, neighbor_dict, to_computed_neighbors,
                        vertex_index_map, edge_weight_map);
        bool change_flag = true;
        while (change_flag) {
            change_flag = false;
            //First Init: Add Neighbor to Check List
            vector<Member> to_check_list;
            for (auto &neighbor_pair:neighbor_dict) { to_check_list.emplace_back(neighbor_pair.second); }
            auto degree_cmp_obj = [this](auto &left, auto &right) {
                return degree(this->vertices_[left.member_index_], *this->graph_ptr_) <
                       degree(this->vertices_[right.member_index_], *this->graph_ptr_);
            };
            sort(to_check_list.begin(), to_check_list.end(), degree_cmp_obj);
            for (auto &neighbor:to_check_list) {
                if (CalDensity(community) < CalDensity(community, neighbor, MutationType::add_neighbor)) {
                    change_flag = true;
                    community.UpdateInfoForMutation(neighbor, MutationType::add_neighbor);
                    neighbor_dict.erase(neighbor.member_index_);
                    auto check_vertex = vertices_[neighbor.member_index_];
                    member_dict.emplace(neighbor.member_index_, neighbor);
                    UpdateStates(check_vertex, community, member_dict, neighbor_dict,
                                 MutationType::add_neighbor, vertex_index_map, edge_weight_map);
                }
            }

            //Second Init: Add Member to Check List
            to_check_list.clear();
            for (auto &member_pair:member_dict) { to_check_list.emplace_back(member_pair.second); }
            sort(to_check_list.begin(), to_check_list.end(), degree_cmp_obj);
            for (auto &member:to_check_list) {
                if (CalDensity(community) < CalDensity(community, member, MutationType::remove_member)) {
                    change_flag = true;
                    community.UpdateInfoForMutation(member, MutationType::add_neighbor);
                    member_dict.erase(member.member_index_);
                    auto check_vertex = vertices_[member.member_index_];
                    neighbor_dict.emplace(member.member_index_, member);
                    UpdateStates(check_vertex, community, member_dict, neighbor_dict,
                                 MutationType::remove_member, vertex_index_map, edge_weight_map);
                }
            }
            community = SplitAndChoose(community.member_indices_);
        }

        auto member_vec = MemberVec(community.member_indices_.size());
        auto local_index = 0;
        for (auto member_index:community.member_indices_) {
            member_vec[local_index] = member_index;
            local_index++;
        }
        //For Later Sort-Merge-Join
        sort(member_vec.begin(), member_vec.end());
        return member_vec;
    }

    double Cis::GetIntersectRatio(MemberVec &left_community, MemberVec &right_community) const {
        auto intersect_set = vector<int>(left_community.size() + right_community.size());
        auto iter_end = set_intersection(left_community.begin(), left_community.end(),
                                         right_community.begin(), right_community.end(), intersect_set.begin());
        auto intersect_set_size = iter_end - intersect_set.begin();
        double rate = static_cast<double>(intersect_set.size()) / min(left_community.size(), right_community.size());
        return rate;
    }

    MemberVec Cis::GetUnion(MemberVec &left_community, MemberVec &right_community) const {
        auto union_set = vector<int>(left_community.size() + right_community.size());
        auto iter_end = set_union(left_community.begin(), left_community.end(),
                                  right_community.begin(), right_community.end(), union_set.begin());
        union_set.resize(iter_end - union_set.begin());
        return union_set;
    }

    void Cis::MergeCommToGlobal(MemberVec &community) {
        if (overlap_community_vec_.size() == 0) {
            overlap_community_vec_.emplace_back(community);
        } else {
            bool is_insert = true;
            for (auto &comm_ptr:overlap_community_vec_) {
                auto cover_rate = GetIntersectRatio(comm_ptr, community);
                if (cover_rate > 1 - DOUBLE_ACCURACY) {
                    comm_ptr = GetUnion(comm_ptr, community);
                    is_insert = false;
                    break;
                }
            }
            if (is_insert) {
                overlap_community_vec_.emplace_back(community);
            }
        }
    }

    Cis::OverlappingCommunityVec Cis::ExecuteCis() {
        auto overlapping_communities_ptr = OverlappingCommunityVec();
        property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
        for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) {
            Vertex vertex = *vp.first;
            auto partial_comm_members = MemberSet();
            partial_comm_members.emplace(vertex_index_map[vertex]);
            auto result_community = ExpandSeed(partial_comm_members);
            MergeCommToGlobal(result_community);
        }
        return overlapping_communities_ptr;
    }
}
