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

#include "util/pretty_print.h"

using namespace std;
using namespace boost;

namespace yche {
    constexpr double DOUBLE_ACCURACY = 0.00001;
    using EntityIdxSet = std::unordered_set<int>;
    using EntityIdxVec = vector<int>;

    struct Entity {
        int entity_index_;
        double w_in_;
        double w_out_;

        Entity(int member_index) : entity_index_(member_index), w_in_(0), w_out_(0) {}
    };

    using EntityDict = std::unordered_map<int, Entity>;

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

        void UpdateInfoForMutation(const Entity &member_info, MutationType mutation_type) {
            if (mutation_type == MutationType::add_neighbor) {
                w_in_ += member_info.w_in_;
                w_out_ += member_info.w_out_;
                member_indices_.emplace(member_info.entity_index_);
            } else {
                w_in_ -= member_info.w_in_;
                w_out_ -= member_info.w_out_;
                member_indices_.emplace(member_info.entity_index_);
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
        using OverlappingCommunityVec=vector<EntityIdxVec>;

        OverlappingCommunityVec overlap_community_vec_;

        Cis(unique_ptr<Graph> &graph_ptr, double lambda);

        OverlappingCommunityVec ExecuteCis();


    private:
        unique_ptr<Graph> graph_ptr_;
        vector<Vertex> vertices_;
        double lambda_;

        double CalDensity(int size, double w_in, double w_out, double lambda) const;

        double CalDensity(Community &community) const;

        double CalDensity(Community &community, Entity &member, MutationType mutation_type) const;

        void InitializeSeeds(const EntityIdxSet &seed, Community &community, EntityDict &member_dict,
                             EntityDict &neighbor_dict, property_map<Graph, vertex_index_t>::type &vertex_index_map,
                             property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateForAddNeighbor(const Vertex &mutate_vertex, Community &community,
                                  EntityDict &member_dict, EntityDict &neighbor_dict,
                                  property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                  property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateForRemoveMember(const Vertex &mutate_vertex, Community &community,
                                   EntityDict &member_dict, EntityDict &neighbor_dict,
                                   property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                   property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void MutateStates(MutationType mutation_type, vector<Entity> to_check_list,
                          Community &community, EntityDict &expand_entity_dict,
                          EntityDict &shrink_entity_dict, auto degree_cmp_obj, bool &change_flag,
                          property_map<Graph, vertex_index_t>::type &vertex_index_map,
                          property_map<Graph, edge_weight_t>::type &edge_weight_map);

        Community FindConnectedComponent(EntityIdxSet &member_set, EntityIdxSet &mark_set, int first_mem_idx,
                                         property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                         property_map<Graph, edge_weight_t>::type &edge_weight_map);

        Community SplitAndChoose(EntityIdxSet &member_set);

        EntityIdxVec ExpandSeed(EntityIdxSet &entity_idx_set);

        double GetIntersectRatio(EntityIdxVec &left_community, EntityIdxVec &right_community) const;

        EntityIdxVec GetUnion(EntityIdxVec &left_community, EntityIdxVec &right_community) const;

        void MergeCommToGlobal(EntityIdxVec &result_community);
    };
}

#endif //CODES_YCHE_CIS_SEQUENTIAL_ALGORITHM_H
