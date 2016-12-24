//
// Created by cheyulin on 5/12/16.
//

#ifndef CODES_YCHE_INPUT_OUTPUT_HANDLER_H
#define CODES_YCHE_INPUT_OUTPUT_HANDLER_H

#include <boost/regex.hpp>
#include <fstream>
#include <sstream>
#include <memory>
#include <map>
#include <vector>

#include "parallel_utils/dataflow_scheduler.h"
#include "parallel_utils/reduce_scheduler.h"


namespace yche {
    using namespace std;

    struct Edge {
        int src_index_;
        int dst_index_;
        double edge_weight_;

        Edge(int src_index_, int dst_index_, double edge_weight_) :
                src_index_(src_index_), dst_index_(dst_index_), edge_weight_(edge_weight_) {}
    };

    vector<Edge> ReadWeightedEdgeList(auto &&file_name_str) {
        auto edges_vec = vector<Edge>();
        auto in_stream = ifstream(file_name_str);
        auto s = string();
        if (!in_stream) {
            cout << "Error opening " << string(file_name_str) << " for input" << endl;
            exit(-1);
        } else {
            while (getline(in_stream, s)) {
                boost::regex pat("#.*edge_weight.*");
                boost::smatch matches;
                cout << s << endl;
                if (boost::regex_match(s, matches, pat))
                    continue;

                auto first_vertex_name = -1;
                auto second_vertex_name = -1;
                auto edge_weight = -1.0;
                auto string_stream = stringstream(s);
                string_stream >> first_vertex_name >> second_vertex_name >> edge_weight;
                edges_vec.emplace_back(first_vertex_name, second_vertex_name, edge_weight);
            }
        }
        return edges_vec;
    }

    vector<pair<int, int>> ReadEdgeList(auto &&file_name_str) {
        auto edges_vec = vector<pair<int, int>>();
        auto in_stream = ifstream(file_name_str);
        auto s = string();
        if (!in_stream) {
            cout << "Error opening " << string(file_name_str) << " for input" << endl;
            exit(-1);
        } else {
            while (getline(in_stream, s)) {
                boost::regex pat("#.*");
                boost::smatch matches;
                cout << s << endl;
                if (boost::regex_match(s, matches, pat))
                    continue;

                auto first_vertex_name = -1;
                auto second_vertex_name = -1;
                auto string_stream = stringstream(s);
                string_stream >> first_vertex_name >> second_vertex_name;
                edges_vec.emplace_back(first_vertex_name, second_vertex_name);
            }
        }
        return edges_vec;
    }

    template<typename Algorithm>
    void ExecuteAlgorithm(int thread_num, unique_ptr<Algorithm> &algorithm_ptr, map<int, int> &index_name_map) {
        cout << "Reduce Enabled" << endl;
        DataFlowScheduler<Algorithm> parallelizer(thread_num, std::move(algorithm_ptr));
        parallelizer.ParallelExecute();
        algorithm_ptr = std::move(parallelizer.algorithm_ptr_);

#ifndef NOT_COUT_COMMUNITY_RESULT
        auto communities_ptr_vec = std::move(algorithm_ptr->overlap_community_vec_);
        cout << "comm_size:" << communities_ptr_vec->size() << endl;
        for (auto &&community_ptr:*communities_ptr_vec) {
            for (auto member_id:*community_ptr) {
                cout << index_name_map[member_id] << "(" << member_id << "),";
            }
            cout << endl;
        }
#endif
    }
}

namespace std {
    ostream &operator<<(ostream &my_ostream, yche::Edge &edge) {
        cout << "src:" << edge.src_index_ << ",dst:" << edge.dst_index_ << ",weight:" << edge.edge_weight_ << endl;
        return my_ostream;
    }
}
#endif //CODES_YCHE_INPUT_OUTPUT_HANDLER_H
