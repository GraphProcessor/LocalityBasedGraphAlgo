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

    void ReadEdgeListWithWeightInToEdgeVector(char *const &file_name_ptr, vector<Edge> &edges_vec) {
        ifstream fin(file_name_ptr);
        string s;
        if (!fin) {
            cout << "Error opening " << string(file_name_ptr) << " for input" << endl;
            exit(-1);
        }

        int i = 0;
        while (getline(fin, s)) {
            using namespace boost;
            boost::regex pat("$#.*edge_weight");
            boost::smatch matches;
            cout << s << endl;
            if (boost::regex_match(s, matches, pat))
                continue;

            int first_vertex_name = -1;
            int second_vertex_name = -1;
            double edge_weight = -1;
            stringstream string_stream;
            string_stream.clear();
            string_stream.str(s);
            string_stream >> first_vertex_name >> second_vertex_name >> edge_weight;
            cout << "src:" << first_vertex_name << ",dst:" << second_vertex_name
                 << ",edge_weight:" << edge_weight << endl;
            edges_vec.emplace_back(first_vertex_name, second_vertex_name, edge_weight);
            i++;
        }
    }

    void ReadEdgeListInToEdgeVector(char *const &file_name_ptr, vector<pair<int, int>> &edges_vec) {
        ifstream fin(file_name_ptr);
        if (!fin) {
            cout << "Error opening " << string(file_name_ptr) << " for input" << endl;
            exit(-1);
        }

        string s;
        int line_num = 0;
        while (getline(fin, s)) {
            using namespace boost;
            boost::regex pat("$#.*");
            boost::smatch matches;
            if (boost::regex_match(s, matches, pat))
                continue;

            auto first_vertex_name = -1;
            auto second_vertex_name = -1;
            stringstream string_stream;
            string_stream.clear();
            string_stream.str(s);
            string_stream >> first_vertex_name;
            string_stream >> second_vertex_name;
            edges_vec.emplace_back(first_vertex_name, second_vertex_name);
            line_num++;
        }
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
#endif //CODES_YCHE_INPUT_OUTPUT_HANDLER_H
