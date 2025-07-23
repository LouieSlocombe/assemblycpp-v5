/**
 * @file vf2.h
 * @brief code to utilise boost vf2 graph isomorphism algorithm
 */
#pragma once
#include <boost/graph/adjacency_list.hpp>            // for adjacency_list
#include <string>                                    // for basic_string
#include <vector>                                    // for vector
#include "boost/graph/graph_selectors.hpp"           // for undirectedS
#include "boost/graph/mcgregor_common_subgraphs.hpp" // for property_map_eq...
#include "boost/graph/properties.hpp"                // for edge_name_t
#include "boost/iterator/iterator_facade.hpp"        // for operator!=
#include "boost/pending/property.hpp"                // for property, no_pr...
#include "globalPrimitives.h"                        // for edgeL, standard...
struct molGraph;
typedef boost::property<boost::edge_name_t, char> bond_vf2;
typedef boost::property<boost::vertex_name_t, std::string, boost::property<boost::vertex_index_t, int>> atom_vf2;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, atom_vf2, bond_vf2> molGraphBoost;
typedef boost::property_map<molGraphBoost, boost::vertex_name_t>::type vertex_name_map_t;
typedef boost::property_map_equivalent<vertex_name_map_t, vertex_name_map_t> vertex_comp_t;
typedef boost::property_map<molGraphBoost, boost::edge_name_t>::type edge_name_map_t;
typedef boost::property_map_equivalent<edge_name_map_t, edge_name_map_t> edge_comp_t;

/**
 * @brief Converts a boolean edgelist based on a standard molGraph and edgeList into a Boost molGraph
 *
 * @param mg The standard molGraph
 * @param edgeList The standard edgeList
 * @param mask The input boolean edgelist
 * @return molGraphBoost the output
 */
molGraphBoost edgelistToBoost(molGraph &mg, std::vector<edgeL> &edgeList, const standardBitset &mask);

/**
 * @brief vf2 graph isomorphism caller
 *
 * @param mmg First graph to be compared
 * @param tmg Second graph to be compared
 * @return true if the two graphs are isomorphic
 * @return false otherwise
 */
bool vf2GraphIso(molGraphBoost &mmg, molGraphBoost &tmg);