//
// Created by cheyulin on 12/20/16.
//

#ifndef CODES_YCHE_CIS_SEQUENTIAL_ALGORITHM_H
#define CODES_YCHE_CIS_SEQUENTIAL_ALGORITHM_H

#include <boost/graph/adjacency_list.hpp>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <vector>
#include <limits>
#include <iostream>

using namespace std;
using namespace boost;

namespace yche {
    constexpr double DOUBLE_ACCURACY = 0.00001;
    using MemberSet = std::unordered_set<int>;
    using MemberVec = std::vector<int>;

    enum class MutationType {
        add_neighbor,
        remove_member
    };

    struct Member {
        int member_index_;
        double w_in_;
        double w_out_;

        Member(int member_index) : member_index_(member_index), w_in_(0), w_out_(0) {}
    };

    using MemberMap = std::unordered_map<int, Member>;

    struct Community {
        std::unordered_set<int> member_indices_;
        double w_in_;
        double w_out_;

        Community(double w_in, double w_out) : w_in_(w_in), w_out_(w_out) {}

        Community(const Community &community) = default;

        Community(Community &&community) : member_indices_(std::move(community.member_indices_)),
                                           w_in_(community.w_in_), w_out_(community.w_out_) {}

        Community &operator=(Community &&community) {
            member_indices_ = std::move(community.member_indices_);
            w_in_ = community.w_in_;
            w_out_ = community.w_out_;
            return *this;
        }

        void UpdateInfoForMutation(const Member &member_info, MutationType mutation_type) {
            if (mutation_type == MutationType::add_neighbor) {
                this->w_in_ += member_info.w_in_;
                this->w_out_ += member_info.w_out_;
                member_indices_.emplace(member_info.member_index_);
            } else {
                this->w_in_ -= member_info.w_in_;
                this->w_out_ -= member_info.w_out_;
                member_indices_.emplace(member_info.member_index_);
            }
        }
    };

    class Cis {
    public:
        using EdgeProperties = property<edge_weight_t, double>;
        using VertexProperties = property<vertex_index_t, int>;
        using Graph = adjacency_list<hash_setS, vecS, undirectedS, VertexProperties, EdgeProperties>;
        using Vertex = graph_traits<Graph>::vertex_descriptor;
        using Edge = graph_traits<Graph>::edge_descriptor;

        using OverlappingCommunityVec=vector<MemberVec>;
        using BasicDataType = MemberSet;
        using MergeDataType = MemberVec;

        OverlappingCommunityVec overlap_community_vec_;

        [[deprecated("Replaced With Parallel Execution")]]
        OverlappingCommunityVec ExecuteCis();

        Cis(unique_ptr<Graph> &graph_ptr, double lambda) : lambda_(lambda), graph_ptr_(std::move(graph_ptr)) {
            property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
            for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) {
                Vertex vertex = *vp.first;
                vertices_.emplace_back(*vp.first);
            }
        }

    private:
        unique_ptr<Graph> graph_ptr_;
        vector<Vertex> vertices_;
        double lambda_;

        double CalDensity(int size, double w_in, double w_out, double lambda) const;

        double CalDensity(Community &community_info_ptr) const;

        double CalDensity(Community &community_info_ptr, Member &member_info_ptr, MutationType mutation_type) const;

        Community SplitAndChoose(MemberSet &member_set);

        void UpdateIterationStates(const Vertex &mutate_vertex,
                                   Community &community_info_ptr,
                                   MemberMap &members,
                                   MemberMap &neighbors, const MutationType &mutation_type,
                                   property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                   property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateForAddNeighbor(const Vertex &mutate_vertex, Community &community_info_ptr,
                                  MemberMap &members, MemberMap &neighbors,
                                  property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                  property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateForRemoveMember(const Vertex &mutate_vertex, Community &community_info_ptr,
                                   MemberMap &members, MemberMap &neighbors,
                                   property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                   property_map<Graph, edge_weight_t>::type &edge_weight_map);


        void InitializeSeeds(const MemberSet &seed_member_ptr,
                             Community &community_info_ptr, MemberMap &members,
                             MemberMap &neighbors, MemberSet &to_computed_neighbors,
                             property_map<Graph, vertex_index_t>::type &vertex_index_map,
                             property_map<Graph, edge_weight_t>::type &edge_weight_map);

        MemberVec ExpandSeed(MemberSet &seed_member_ptr);

        double GetIntersectRatio(MemberVec &left_community, MemberVec &right_community);

        MemberVec MergeBoth(MemberVec &left_community, MemberVec &right_community);

        void MergeCommToGlobal(OverlappingCommunityVec &community_collection, MergeDataType &result);
    };
}


#endif //CODES_YCHE_CIS_SEQUENTIAL_ALGORITHM_H
