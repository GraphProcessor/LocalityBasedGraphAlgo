//
// Created by cheyulin on 4/15/16.
//

#ifndef CODES_YCHE_CIS_ALGORITHM_H
#define CODES_YCHE_CIS_ALGORITHM_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/range.hpp>
#include <memory>
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

    using IndexType = unsigned long;
    using CommunityMemberSet = std::unordered_set<IndexType>;
    using CommunityMemberVec = std::vector<IndexType>;

    enum class MutationType {
        add_neighbor,
        remove_member
    };


    //Cache the Previous Computation Results for a Single Vertex
    struct MemberInfo {
        IndexType member_index_;
        double w_in_;
        double w_out_;

        MemberInfo(IndexType member_index) : member_index_(member_index), w_in_(0), w_out_(0) {}
    };

    using MemberInfoMap = std::unordered_map<IndexType, unique_ptr<MemberInfo>>;

    //Cache the Previous Computation Results for a Community
    struct CommunityInfo {
        unique_ptr<CommunityMemberSet> members_;
        double w_in_;
        double w_out_;

        CommunityInfo(double w_in, double w_out) : w_in_(w_in), w_out_(w_out) {}

        void UpdateInfoForMutation(const MemberInfo &member_info, MutationType mutation_type) {
            if (mutation_type == MutationType::add_neighbor) {
                this->w_in_ += member_info.w_in_;
                this->w_out_ += member_info.w_out_;
                this->members_->insert(member_info.member_index_);
            } else {
                this->w_in_ -= member_info.w_in_;
                this->w_out_ -= member_info.w_out_;
                this->members_->erase(member_info.member_index_);
            }
        }
    };

    class Cis {
    public:
        using EdgeProperties = property<edge_weight_t, double>;
        using VertexProperties = property<vertex_index_t, IndexType>;
        using Graph = adjacency_list<hash_setS, vecS, undirectedS, VertexProperties, EdgeProperties>;
        using Vertex = graph_traits<Graph>::vertex_descriptor;
        using Edge = graph_traits<Graph>::edge_descriptor;

        using OverlappingCommunityVec=vector<unique_ptr<CommunityMemberVec>>;

        //Start Implementation Interfaces For DataFlower Traits
        unique_ptr<OverlappingCommunityVec> overlap_community_vec_;
        using BasicDataType = CommunityMemberSet;
        using MergeDataType = CommunityMemberVec;

        unique_ptr<vector<unique_ptr<BasicDataType>>> InitBasicComputationData();

        unique_ptr<MergeDataType> LocalComputation(unique_ptr<BasicDataType> seed_member_ptr);

        void MergeToGlobal(unique_ptr<MergeDataType> &result);

        //Start Implementation Interfaces For Reducer Traits
        using ReduceDataType = OverlappingCommunityVec;

        unique_ptr<ReduceDataType> WrapMergeDataToReduceData(unique_ptr<MergeDataType> &merge_data_ptr);

        function<bool(unique_ptr<ReduceDataType> &, unique_ptr<ReduceDataType> &)> CmpReduceData;

        function<unique_ptr<ReduceDataType>(unique_ptr<ReduceDataType> &,
                                            unique_ptr<ReduceDataType> &right_data_ptr)> ReduceComputation;

        //Start Implementation Interfaces For Fine-Grained-Merge-Scheduler Traits
        using ElementReferenceType = typename boost::range_reference<ReduceDataType>::type;

        function<bool(ElementReferenceType, ElementReferenceType)> PairMergeComputation;

        function<void(ElementReferenceType, ElementReferenceType)> SuccessAction;

        function<void(ElementReferenceType, unique_ptr<ReduceDataType> &)> FailAction;

        [[deprecated("Replaced With Parallel Execution")]]
        unique_ptr<OverlappingCommunityVec> ExecuteCis();

        Cis(unique_ptr<Graph> graph_ptr, double lambda) : lambda_(lambda), graph_ptr_(std::move(graph_ptr)) {
            vertices_.clear();
            //Init Vertices
            property_map<Graph, vertex_index_t>::type vertex_index_map = boost::get(vertex_index, *graph_ptr_);
            for (auto vp = vertices(*graph_ptr_); vp.first != vp.second; ++vp.first) {
                Vertex vertex = *vp.first;
                vertices_.push_back(*vp.first);
            }

            overlap_community_vec_ = make_unique<OverlappingCommunityVec>();

            CmpReduceData = [](unique_ptr<ReduceDataType> &left, unique_ptr<ReduceDataType> &right) -> bool {
                auto cmp = [](auto &tmp_left, auto &tmp_right) -> bool {
                    return tmp_left->size() < tmp_right->size();
                };
                auto iter1 = max_element(left->begin(), left->end(), cmp);
                auto iter2 = max_element(left->begin(), left->end(), cmp);
                return (*iter1)->size() > (*iter2)->size();
            };

            ReduceComputation = [this](unique_ptr<ReduceDataType> &left_data_ptr,
                                       unique_ptr<ReduceDataType> &right_data_ptr) -> unique_ptr<ReduceDataType> {
                for (auto &right_merge_data:*right_data_ptr) {
                    MergeToCommunityCollection(left_data_ptr, right_merge_data);
                }
                return std::move(left_data_ptr);
            };

            //Start Implementation Interfaces For Fine-Grained-Merge-Scheduler Traits
            PairMergeComputation = [this](ElementReferenceType left_element_ptr,
                                          ElementReferenceType right_element_ptr) -> bool {
                return GetTwoCommunitiesCoverRate(left_element_ptr, right_element_ptr) > 1 - DOUBLE_ACCURACY;
            };

            SuccessAction = [this](ElementReferenceType left_element_ptr, ElementReferenceType right_element_ptr) {
                right_element_ptr = std::move(MergeTwoCommunities(right_element_ptr, left_element_ptr));
            };

            FailAction = [](ElementReferenceType left_element_ptr, unique_ptr<ReduceDataType> &reduce_data_ptr) {
                reduce_data_ptr->push_back(std::move(left_element_ptr));
            };
        }

    private:

        unique_ptr<Graph> graph_ptr_;
        vector<Vertex> vertices_;

        double lambda_;

        double CalculateDensity(const IndexType &size, const double &w_in, const double &w_out, const double &lambda);

        inline double CalculateDensity(unique_ptr<CommunityInfo> &community_info_ptr);

        inline double CalculateDensity(unique_ptr<CommunityInfo> &community_info_ptr,
                                       unique_ptr<MemberInfo> &member_info_ptr, const MutationType &mutation_type);

        unique_ptr<CommunityInfo> SplitAndChooseBestConnectedComponent(unique_ptr<CommunityMemberSet> &community_ptr);

        void UpdateMembersNeighborsCommunityInfo(const Vertex &mutate_vertex,
                                                 unique_ptr<CommunityInfo> &community_info_ptr,
                                                 MemberInfoMap &members,
                                                 MemberInfoMap &neighbors, const MutationType &mutation_type,
                                                 property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                                 property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateMembersNeighborsCommunityInfoForAddNeighbor(const Vertex &mutate_vertex,
                                                               unique_ptr<CommunityInfo> &community_info_ptr,
                                                               MemberInfoMap &members,
                                                               MemberInfoMap &neighbors,
                                                               property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                                               property_map<Graph, edge_weight_t>::type &edge_weight_map);

        void UpdateMembersNeighborsCommunityInfoForRemoveMember(const Vertex &mutate_vertex,
                                                                unique_ptr<CommunityInfo> &community_info_ptr,
                                                                MemberInfoMap &members,
                                                                MemberInfoMap &neighbors,
                                                                property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                                                property_map<Graph, edge_weight_t>::type &edge_weight_map);


        void InitializationForSeedExpansion(const unique_ptr<CommunityMemberSet> &seed_member_ptr,
                                            unique_ptr<CommunityInfo> &community_info_ptr, MemberInfoMap &members,
                                            MemberInfoMap &neighbors, CommunityMemberSet &to_computed_neighbors,
                                            property_map<Graph, vertex_index_t>::type &vertex_index_map,
                                            property_map<Graph, edge_weight_t>::type &edge_weight_map);

        unique_ptr<CommunityMemberVec> ExpandSeed(unique_ptr<CommunityMemberSet> &seed_member_ptr);

        double GetTwoCommunitiesCoverRate(unique_ptr<CommunityMemberVec> &left_community,
                                          unique_ptr<CommunityMemberVec> &right_community);

        unique_ptr<CommunityMemberVec> MergeTwoCommunities(unique_ptr<CommunityMemberVec> &left_community,
                                                           unique_ptr<CommunityMemberVec> &right_community);

        void MergeToCommunityCollection(decltype(overlap_community_vec_) &community_collection,
                                        unique_ptr<MergeDataType> &result);
    };
}

#endif //CODES_YCHE_CIS_ALGORITHM_H
