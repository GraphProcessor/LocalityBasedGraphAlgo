//
// Created by cheyulin on 12/20/16.
//

#include "cis_sequential_algorithm.h"


namespace yche {
    double Cis::CalDensity(int size, double w_in, double w_out, double lambda) const {
        if (size < 1)
            return numeric_limits<double>::min();
        else {
            double partA = ((1 - lambda) * (w_in / (w_in + w_out)));
            double partB = (lambda * ((2 * w_in) / (size * (size - 1))));
            if (size == 1)
                partB = lambda;
            return partA + partB;
        }
    }

    double Cis::CalDensity(Community &community_info_ptr) const {
        return CalDensity(static_cast<int>(community_info_ptr.member_indices_.size()),
                          community_info_ptr.w_in_,
                          community_info_ptr.w_out_, lambda_);
    }

    double Cis::CalDensity(Community &community_info_ptr,
                           Member &member_info_ptr, MutationType mutation_type) const {
        if (mutation_type == MutationType::add_neighbor)
            return CalDensity(static_cast<int>(community_info_ptr.member_indices_.size() + 1),
                              community_info_ptr.w_in_ + member_info_ptr.w_in_,
                              community_info_ptr.w_out_ + member_info_ptr.w_out_, lambda_);
        else
            return CalDensity(static_cast<int>(community_info_ptr.member_indices_.size() - 1),
                              community_info_ptr.w_in_ - member_info_ptr.w_in_,
                              community_info_ptr.w_out_ - member_info_ptr.w_out_, lambda_);

    }

    Community Cis::SplitAndChoose(MemberSet &member_set) {
        queue<int> frontier;
        std::unordered_set<int> mark_set;
        vector<Community> community_info_ptr_vec;
        while (member_set.size() > 0) {
            int first_vertex = *member_set.begin();
            //Enqueue
            frontier.push(first_vertex);
            mark_set.emplace(first_vertex);

            auto community_info_ptr = Community(0, 0);
            //One Connected Component
            while (frontier.size() > 0) {
                property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
                property_map<Graph, edge_weight_t>::type edge_weight_map = boost::get(edge_weight, *graph_ptr_);
                Edge edge;
                bool edge_exist_flag;
                auto expand_vertex_index = frontier.front();
                auto expand_vertex = vertices_[expand_vertex_index];

                //Add To Community Only At Start of Expansion
                community_info_ptr.member_indices_.emplace(expand_vertex_index);
                for (auto vp = adjacent_vertices(vertices_[expand_vertex_index], *graph_ptr_);
                     vp.first != vp.second; ++vp.first) {
                    auto neighbor_vertex = *vp.first;
                    auto adjacency_vertex_index = static_cast<int>(vertex_index_map[neighbor_vertex]);
                    tie(edge, edge_exist_flag) = boost::edge(expand_vertex, neighbor_vertex, *graph_ptr_);
                    //Do W_in and W_out Computation
                    auto iter = member_set.find(adjacency_vertex_index);
                    if (mark_set.find(adjacency_vertex_index) == mark_set.end() && iter != member_set.end()) {
                        community_info_ptr.w_in_ += edge_weight_map[edge];
                        //Erase and Enqueue
                        frontier.push(adjacency_vertex_index);
                        mark_set.emplace(adjacency_vertex_index);
                    } else {
                        community_info_ptr.w_out_ = edge_weight_map[edge];
                    }
                }
                member_set.erase(expand_vertex_index);
                frontier.pop();
            }
            community_info_ptr_vec.emplace_back(std::move(community_info_ptr));
        }

        //Sort and Select Best Community , i.e., With Largest Density
        sort(community_info_ptr_vec.begin(), community_info_ptr_vec.end(),
             [this](auto &left_comm_info_ptr, auto &right_comm_info_ptr) {
                 double left_density = this->CalDensity(left_comm_info_ptr);
                 double right_density = this->CalDensity(right_comm_info_ptr);
                 if (left_density != right_density) {
                     return left_density > right_density;
                 } else {
                     return (left_comm_info_ptr.member_indices_).size() > (right_comm_info_ptr.member_indices_).size();
                 }
             });
        return community_info_ptr_vec[0];
    }

    void Cis::InitializeSeeds(const MemberSet &seed_member_ptr, Community &community_info_ptr,
                              MemberMap &members, MemberMap &neighbors,
                              MemberSet &to_computed_neighbors,
                              property_map<Graph, vertex_index_t>::type &vertex_index_map,
                              property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        //Tally Members of the seed, calculating individual w_in and w_out
        for (auto &seed_vertex_index :seed_member_ptr) {
            auto member_info_ptr = Member(seed_vertex_index);
            community_info_ptr.member_indices_.emplace(seed_vertex_index);
            Vertex seed_vertex = vertices_[seed_vertex_index];

            for (auto vp = adjacent_vertices(seed_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                auto neighbor_vertex_index = static_cast<int>(vertex_index_map[*vp.first]);
                auto edge_weight = edge_weight_map[
                        edge(seed_vertex, vertices_[neighbor_vertex_index], *graph_ptr_).first];
                if (seed_member_ptr.find(neighbor_vertex_index) != seed_member_ptr.end()) {
                    member_info_ptr.w_in_ += edge_weight;
                    community_info_ptr.w_in_ += edge_weight;
                } else {
                    member_info_ptr.w_out_ += edge_weight;
                    community_info_ptr.w_out_ += edge_weight;
                    to_computed_neighbors.emplace(neighbor_vertex_index);
                }
            }
            members.emplace(member_info_ptr.member_index_, member_info_ptr);
        }

        //Tally For Neighbors of the seed, calculate w_int and w_out
        for (auto &neighbor_vertex_index :to_computed_neighbors) {
            auto neighbor_info_ptr = Member(neighbor_vertex_index);
            Vertex neighbor_vertex = vertices_[neighbor_vertex_index];
            for (auto vp = adjacent_vertices(neighbor_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
                auto neighbor_neighbor_vertex_index = static_cast<int>(vertex_index_map[*vp.first]);
                auto edge_weight = edge_weight_map[edge(neighbor_vertex, vertices_[neighbor_neighbor_vertex_index],
                                                        *graph_ptr_).first];
                if (seed_member_ptr.find(neighbor_neighbor_vertex_index) != seed_member_ptr.end()) {
                    neighbor_info_ptr.w_in_ += edge_weight;
                } else {
                    neighbor_info_ptr.w_out_ += edge_weight;
                }
            }
            neighbors.emplace(neighbor_info_ptr.member_index_, neighbor_info_ptr);
        }
    }

    void Cis::UpdateIterationStates(const Cis::Vertex &mutate_vertex,
                                    Community &community_info_ptr, MemberMap &members, MemberMap &neighbors,
                                    const MutationType &mutation_type,
                                    property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                    property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        if (mutation_type == MutationType::add_neighbor) {
            UpdateForAddNeighbor(mutate_vertex, community_info_ptr, members,
                                 neighbors, vertex_index_map, edge_weight_map);
        } else {
            UpdateForRemoveMember(mutate_vertex, community_info_ptr, members,
                                  neighbors, vertex_index_map, edge_weight_map);
        }
    }

    void Cis::UpdateForAddNeighbor(const Cis::Vertex &mutate_vertex,
                                   Community &community_info_ptr,
                                   MemberMap &members, MemberMap &neighbors,
                                   property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                   property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        //Update Member and Neighbor List
        for (auto vp = adjacent_vertices(mutate_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto check_neighbor_vertex = *vp.first;
            auto check_neighbor_vertex_index = vertex_index_map[check_neighbor_vertex];
            auto check_neighbor_ptr = make_unique<Member>(check_neighbor_vertex_index);
            auto edge_weight = edge_weight_map[edge(mutate_vertex, check_neighbor_vertex,
                                                    *graph_ptr_).first];

            auto iter = members.find(check_neighbor_ptr->member_index_);
            if (iter != members.end() ||
                (iter = neighbors.find(check_neighbor_ptr->member_index_)) != neighbors.end()) {
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
                                                       vertices_[neighbor_neighbor_vertex_index],
                                                       *graph_ptr_).first];
                    if (community_info_ptr.member_indices_.find(neighbor_neighbor_vertex_index) !=
                        community_info_ptr.member_indices_.end()) {
                        member_info_ptr.w_in_ += edge_weight;
                    } else {
                        member_info_ptr.w_out_ += edge_weight;
                    }
                }
                neighbors.emplace(member_info_ptr.member_index_, member_info_ptr);
            }
        }
    }

    void Cis::UpdateForRemoveMember(const Cis::Vertex &mutate_vertex,
                                    Community &community_info_ptr,
                                    MemberMap &members, MemberMap &neighbors,
                                    property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                    property_map<Graph, edge_weight_t>::type &edge_weight_map) {
        //Update Member and Neighbor List
        for (auto vp = adjacent_vertices(mutate_vertex, *graph_ptr_); vp.first != vp.second; ++vp.first) {
            auto check_neighbor_vertex = *vp.first;
            auto check_neighbor_vertex_index = vertex_index_map[check_neighbor_vertex];
            auto check_neighbor_ptr = std::move(make_unique<Member>(check_neighbor_vertex_index));
            auto edge_weight = edge_weight_map[edge(mutate_vertex, check_neighbor_vertex,
                                                    *graph_ptr_).first];

            auto iter = members.find(check_neighbor_ptr->member_index_);
            if (iter != members.end() ||
                (iter = neighbors.find(check_neighbor_ptr->member_index_)) != neighbors.end()) {
                //Update Info In Members and Neighbors
                iter->second.w_in_ -= edge_weight;
                iter->second.w_out_ += edge_weight;
            }
        }

    }

    MemberVec Cis::ExpandSeed(MemberSet &seed_member_ptr) {
        //First: InitializeSeeds members and neighbors(or frontier)
        auto community_info_ptr = Community(0, 0);
        MemberMap members;
        MemberMap neighbors;
        MemberSet to_computed_neighbors;
        property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
        property_map<Graph, edge_weight_t>::type edge_weight_map = boost::get(edge_weight, *graph_ptr_);

        InitializeSeeds(seed_member_ptr, community_info_ptr, members, neighbors, to_computed_neighbors,
                        vertex_index_map, edge_weight_map);
        bool change_flag = true;
        while (change_flag) {
            change_flag = false;
            //First Init: Add Neighbor to Check List
            vector<Member> to_check_list;
            for (auto &neighbor_info_ptr:neighbors) { to_check_list.emplace_back(neighbor_info_ptr.second); }
            auto degree_cmp = [this](auto &left_ptr, auto &right_ptr) -> bool {
                return degree(this->vertices_[left_ptr.member_index_], *this->graph_ptr_) <
                       degree(this->vertices_[right_ptr.member_index_], *this->graph_ptr_);
            };
            sort(to_check_list.begin(), to_check_list.end(), degree_cmp);
            //First For Add-Neighbor Iteration, Check all the neighbors in the frontier
            for (auto &neighbor_info_ptr:to_check_list) {
                if (CalDensity(community_info_ptr)
                    < CalDensity(community_info_ptr, neighbor_info_ptr, MutationType::add_neighbor)) {
                    //Add Neighbor
                    change_flag = true;
                    community_info_ptr.UpdateInfoForMutation(neighbor_info_ptr, MutationType::add_neighbor);
                    neighbors.erase(neighbor_info_ptr.member_index_);
                    auto check_vertex = vertices_[neighbor_info_ptr.member_index_];
                    members.emplace(neighbor_info_ptr.member_index_, neighbor_info_ptr);
                    //Update Member and Neighbor List
                    UpdateIterationStates(check_vertex, community_info_ptr, members,
                                          neighbors, MutationType::add_neighbor,
                                          vertex_index_map, edge_weight_map);
                }
            }

            //Second Init: Add Member to Check List
            to_check_list.clear();
            for (auto &member_info_ptr:members) { to_check_list.emplace_back(member_info_ptr.second); }
            sort(to_check_list.begin(), to_check_list.end(), degree_cmp);
            //Second For Remove-Member Iteration, check all the members in the current community
            for (auto &member_info_ptr:to_check_list) {
                if (CalDensity(community_info_ptr)
                    < CalDensity(community_info_ptr, member_info_ptr, MutationType::remove_member)) {
                    //Remove Member
                    change_flag = true;
                    community_info_ptr.UpdateInfoForMutation(member_info_ptr, MutationType::add_neighbor);
                    members.erase(member_info_ptr.member_index_);
                    auto check_vertex = vertices_[member_info_ptr.member_index_];
                    neighbors.emplace(member_info_ptr.member_index_, std::move(member_info_ptr));
                    //Update Member and Neighbor List
                    UpdateIterationStates(check_vertex, community_info_ptr, members,
                                          neighbors, MutationType::remove_member, vertex_index_map, edge_weight_map);
                }
            }
//            community_info_ptr = SplitAndChoose(community_info_ptr.member_indices_);
        }
        MemberVec community_member_vector;
        community_member_vector.resize(community_info_ptr.member_indices_.size());
        auto local_index = 0;
        for (auto member_index:community_info_ptr.member_indices_) {
            community_member_vector[local_index] = member_index;
            local_index++;
        }
        //For Later Sort-Merge-Join
        sort(community_member_vector.begin(), community_member_vector.end());
        return community_member_vector;
    }

    double Cis::GetIntersectRatio(MemberVec &left_community,
                                  MemberVec &right_community) {
        vector<int> intersect_set(left_community.size() + right_community.size());
        auto iter_end = set_intersection(left_community.begin(), left_community.end(),
                                         right_community.begin(), right_community.end(), intersect_set.begin());
        auto intersect_set_size = iter_end - intersect_set.begin();
        double rate = static_cast<double>(intersect_set.size()) / min(left_community.size(), right_community.size());
        return rate;
    }

    MemberVec Cis::MergeBoth(MemberVec &left_community, MemberVec &right_community) {
        vector<int> union_set(left_community.size() + right_community.size());
        auto iter_end = set_union(left_community.begin(), left_community.end(),
                                  right_community.begin(), right_community.end(), union_set.begin());
        union_set.resize(iter_end - union_set.begin());
        return union_set;
    }

    void Cis::MergeCommToGlobal(OverlappingCommunityVec &community_collection, MergeDataType &result) {
        if (community_collection.size() == 0) {
            community_collection.emplace_back(result);
        } else {
            bool is_insert = true;
            for (auto &comm_ptr:community_collection) {
                auto cover_rate = GetIntersectRatio(comm_ptr, result);
                if (cover_rate > 1 - DOUBLE_ACCURACY) {
                    comm_ptr = MergeBoth(comm_ptr, result);
                    is_insert = false;
                    break;
                }
            }
            if (is_insert) {
                community_collection.emplace_back(result);
            }
        }
    }

    [[deprecated("Replaced With Parallel Execution")]]
    Cis::OverlappingCommunityVec Cis::ExecuteCis() {
        auto overlapping_communities_ptr = OverlappingCommunityVec();
        property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
        for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) {
            Vertex vertex = *vp.first;
            auto partial_comm_members = MemberSet();
            partial_comm_members.emplace(vertex_index_map[vertex]);
            auto result_community = ExpandSeed(partial_comm_members);
            MergeCommToGlobal(overlap_community_vec_, result_community);
        }
        return overlapping_communities_ptr;
    }
}
