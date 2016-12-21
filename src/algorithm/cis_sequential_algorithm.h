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
    using MemberIdxSet = std::unordered_set<int>;
    using MemberIdxVec = std::vector<int>;

    struct Member {
        int member_index_;
        double w_in_;
        double w_out_;

        Member(int member_index) : member_index_(member_index), w_in_(0), w_out_(0) {}
    };

    using MemberDict = std::unordered_map<int, Member>;

    enum class MutationType {
        add_neighbor,
        remove_member
    };

    struct Community {
        std::unordered_set<int> member_indices_;
        double w_in_;
        double w_out_;

        Community() : w_in_(0), w_out_(0) {}

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
        using OverlappingCommunityVec=vector<MemberIdxVec>;

        OverlappingCommunityVec overlap_community_vec_;

        Cis(unique_ptr<Graph> &graph_ptr, double lambda);

        OverlappingCommunityVec ExecuteCis();


    private:
        unique_ptr<Graph> graph_ptr_;
        vector<Vertex> vertices_;
        double lambda_;

        double CalDensity(int size, double w_in, double w_out, double lambda) const;

        double CalDensity(Community &community) const;

        double CalDensity(Community &community, Member &member, MutationType mutation_type) const;

        void InitializeSeeds(const MemberIdxSet &seed, Community &community, MemberDict &member_dict,
                             MemberDict &neighbor_dict, property_map<Graph, vertex_index_t>::type &vertex_index_map,
                             property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateForAddNeighbor(const Vertex &mutate_vertex, Community &community,
                                  MemberDict &member_dict, MemberDict &neighbor_dict,
                                  property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                  property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateForRemoveMember(const Vertex &mutate_vertex, Community &community,
                                   MemberDict &member_dict, MemberDict &neighbor_dict,
                                   property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                   property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void MutateStates(MutationType mutation_type, vector<Member> to_check_list,
                          Community &community, MemberDict &expand_member_dict,
                          MemberDict &shrink_member_dict, auto degree_cmp_obj, bool &change_flag,
                          property_map<Graph, vertex_index_t>::type &vertex_index_map,
                          property_map<Graph, edge_weight_t>::type &edge_weight_map);

        Community SplitAndChoose(MemberIdxSet &member_set);

        MemberIdxVec ExpandSeed(MemberIdxSet &seed);

        double GetIntersectRatio(MemberIdxVec &left_community, MemberIdxVec &right_community) const;

        MemberIdxVec GetUnion(MemberIdxVec &left_community, MemberIdxVec &right_community) const;

        void MergeCommToGlobal(MemberIdxVec &result_community);
    };
}

#endif //CODES_YCHE_CIS_SEQUENTIAL_ALGORITHM_H
