//
// Created by cheyulin on 5/16/16.
//

#ifndef OCD_EVALUATION_YCHE_MODULARITY_H
#define OCD_EVALUATION_YCHE_MODULARITY_H

#include "input_output_handler.h"

namespace yche {
    using namespace boost;
    using namespace std;
    constexpr double PRECISION = 0.0001;

    template<typename CoefficientFunc>
    class ModularityLinkBelonging {
    private:
        using Graph = adjacency_list<setS, vecS, undirectedS>;
        using Vertex = graph_traits<Graph>::vertex_descriptor;
        using VertexIndexType = unsigned long;

        unique_ptr<Graph> graph_ptr_;
        map<VertexIndexType, Vertex> name_vertex_map_;
        vector<pair<VertexIndexType, VertexIndexType>> edge_vec_;
        vector<unique_ptr<vector<Vertex>>> overlap_comms_;

        CoefficientFunc coefficient_calc_func_;
        vector<int> community_belongings_count_vec_;
        vector<double> comm_belong_count_rev_vec_;
        vector<VertexIndexType> in_degree_vec_;
        vector<VertexIndexType> out_degree_vec_;

        void ConstructGraph(char *input_file_str);

        void InitAlphaVertexCommunity();

    public:
        ModularityLinkBelonging(char *input_file_str, char *output_file_str, CoefficientFunc coefficient_calc_func_) :
                coefficient_calc_func_(coefficient_calc_func_) {
            ConstructGraph(input_file_str);
            yche::GetCommunities<VertexIndexType, Vertex>(output_file_str, overlap_comms_, name_vertex_map_);
            InitAlphaVertexCommunity();
        }

        double CalculateModularity();
    };

    template<typename CoefficientFunc>
    void ModularityLinkBelonging<CoefficientFunc>::ConstructGraph(char *input_file_str) {
        graph_ptr_ = make_unique<Graph>();
        yche::ReadEdgeList<VertexIndexType>(input_file_str, edge_vec_);

        for (auto &edge:edge_vec_) {
            if (name_vertex_map_.find(edge.first) == name_vertex_map_.end()) {
                Vertex vertex = add_vertex(*graph_ptr_);
                name_vertex_map_.emplace(edge.first, vertex);
            }
            if (name_vertex_map_.find(edge.second) == name_vertex_map_.end()) {
                Vertex vertex = add_vertex(*graph_ptr_);
                name_vertex_map_.emplace(edge.second, vertex);
            }
            add_edge(name_vertex_map_[edge.first], name_vertex_map_[edge.second], *graph_ptr_);
        }
        community_belongings_count_vec_.resize(num_vertices(*graph_ptr_), 0);
        comm_belong_count_rev_vec_.resize(community_belongings_count_vec_.size(), 0);
    }

    template<typename CoefficientFunc>
    double ModularityLinkBelonging<CoefficientFunc>::CalculateModularity() {
        auto modularity_rate = 0.0;
        auto edge_num = num_edges(*graph_ptr_);
        auto vertices_num = num_vertices(*graph_ptr_);
        auto vertex_index_map = get(vertex_index, *graph_ptr_);

        cout << edge_num << "," << vertices_num << endl;
        cout << overlap_comms_.size() << endl;

        for (auto &community_ptr :overlap_comms_) {
            auto &community_members = *community_ptr;
            auto comm_size = community_ptr->size();
            auto f_value_matrix = vector<vector<double>>(comm_size);
            for (auto &row:f_value_matrix) {
                row.resize(comm_size, 0);
            }
            auto f_sum_in_arr = vector<double>(comm_size, 0);
            auto f_sum_out_arr = vector<double>(comm_size, 0);
            auto in_degree_arr = vector<double>(comm_size, 0);
            auto out_degree_arr = vector<double>(comm_size, 0);

            auto adjacent_flag = false;
            //Calculate F_value_matrix F_sum_in F_sum_out
            for (auto i = 0; i < comm_size; i++) {
                auto start_index = vertex_index_map[community_members[i]];
                in_degree_arr[i] = in_degree_vec_[start_index];
                out_degree_arr[i] = out_degree_vec_[start_index];
                for (auto j = 0; j < comm_size; j++) {
                    auto end_index = vertex_index_map[community_members[j]];
                    if (i != j) {
                        adjacent_flag = edge(community_members[i], community_members[j], *graph_ptr_).second;
                        if (adjacent_flag) {
                            f_value_matrix[i][j] = coefficient_calc_func_(
                                    comm_belong_count_rev_vec_[start_index], comm_belong_count_rev_vec_[end_index]);
                            f_sum_out_arr[i] += f_value_matrix[i][j];
                            f_sum_in_arr[j] += f_value_matrix[i][j];
                        }
                    }
                }
            }
            for (auto i = 0; i < comm_size; i++) {
                f_sum_in_arr[i] /= vertices_num;
                f_sum_out_arr[i] /= vertices_num;
            }
            for (auto i = 0; i < comm_size; i++) {
                for (auto j = 0; j < comm_size; j++) {
                    if (i != j && f_value_matrix[i][j] > PRECISION) {
                        modularity_rate += f_value_matrix[i][j] -
                                           out_degree_arr[i] * in_degree_arr[j] * f_sum_out_arr[i] * f_sum_in_arr[j] /
                                           edge_num;
                    }
                }
            }
            cout << "Curr Modularity:" << modularity_rate / edge_num << endl << endl;
        }
        modularity_rate = modularity_rate / edge_num;
        return modularity_rate;
    }

    template<typename CoefficientFunc>
    void ModularityLinkBelonging<CoefficientFunc>::InitAlphaVertexCommunity() {
        auto vertex_index_map = get(vertex_index, *graph_ptr_);
        long long sum = 0;
        for (auto &community_ptr :overlap_comms_) {
            for (auto &vertex :*community_ptr) {
                community_belongings_count_vec_[vertex_index_map[vertex]]++;
            }
        }

        for (auto i = 0; i < community_belongings_count_vec_.size(); i++) {
            if (community_belongings_count_vec_[i] != 0) {
                comm_belong_count_rev_vec_[i] = 1.0 / community_belongings_count_vec_[i];
                sum += community_belongings_count_vec_[i];
            }
        }

        cout << "sum" << sum << endl;
        auto num_vertices = boost::num_vertices(*graph_ptr_);
        in_degree_vec_.resize(num_vertices);
        out_degree_vec_.resize(num_vertices);
        auto graph = *graph_ptr_;

        for (auto iter = vertices(graph); iter.first != iter.second; iter.first++) {
            auto vertex = *iter.first;
            in_degree_vec_[vertex_index_map[vertex]] = in_degree(vertex, graph);
            out_degree_vec_[vertex_index_map[vertex]] = out_degree(vertex, graph);
        }
    }
}

#endif //OCD_EVALUATION_YCHE_MODULARITY_H
