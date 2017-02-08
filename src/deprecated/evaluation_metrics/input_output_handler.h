//
// Created by cheyulin on 5/16/16.
//

#ifndef OCD_EVALUATION_YCHE_INPUT_OUTPUT_HANDLER_H
#define OCD_EVALUATION_YCHE_INPUT_OUTPUT_HANDLER_H

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <cmath>

#include <boost/graph/adjacency_list.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

namespace yche {
    using namespace std;
    using namespace boost;
    struct LinkBelongingModularity;
    struct NodeBelongingModularity;

    template<typename VertexIndexType>
    void ReadEdgeList(const char *file_name_ptr, vector<pair<VertexIndexType, VertexIndexType>> &edges_vec) {
        auto fin = ifstream(file_name_ptr);
        auto s = string();
        if (!fin) {
            cout << "Error opening " << string(file_name_ptr) << " for input" << endl;
            exit(-1);
        } else {
            while (getline(fin, s)) {
                using namespace boost;
                boost::regex pat("$#.*");
                boost::smatch matches;
                if (boost::regex_match(s, matches, pat))
                    continue;

                auto first_vertex_name = -1;
                auto second_vertex_name = -1;
                auto string_stream = stringstream(s);
                string_stream >> first_vertex_name >> second_vertex_name;
                edges_vec.emplace_back(first_vertex_name, second_vertex_name);
            }
        }
    }

    template<typename VertexIndexType, typename VertexType>
    void GetCommunities(const char *output_file, vector<unique_ptr<vector<VertexIndexType>>> &overlap_communities,
                        map<VertexIndexType, VertexType> &name_index_map) {
        ifstream fin(output_file);
        string s;
        if (!fin) {
            cout << "Error opening " << string(output_file) << " for input" << endl;
            exit(-1);
        } else {
            boost::regex pat_start("comm_size.*");
            boost::regex pat_num_list("[0-9]+.*");
            boost::smatch matches;
            bool is_begin_read = false;
            while (getline(fin, s)) {
                using namespace boost;
                if (!is_begin_read) {
                    if (boost::regex_match(s, matches, pat_start))
                        is_begin_read = true;
                    else
                        continue;
                } else if (boost::regex_match(s, matches, pat_num_list)) {
                    boost::regex re(",");
                    s = boost::regex_replace(s, re, " ");
                    auto str_stream = stringstream(s);
                    auto community_ptr = make_unique<vector<VertexIndexType>>();
                    VertexIndexType tmp;
                    while (!str_stream.eof()) {
                        str_stream >> tmp;
                        community_ptr->emplace_back(name_index_map[tmp]);
                    }
                    overlap_communities.push_back(std::move(community_ptr));
                }
            }
        }
    }
}

#endif //OCD_EVALUATION_YCHE_INPUT_OUTPUT_HANDLER_H
