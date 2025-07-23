#include "vf2.h"
#include <algorithm>                             // for copy
#include <boost/graph/adjacency_list.hpp>        // for target, source
#include <boost/graph/vf2_sub_graph_iso.hpp>     // for vertex_order_by_mult
#include <cstddef>                               // for size_t, std
#include <list>                                  // for operator==
#include <map>                                   // for operator==
#include <string>                                // for string
#include <unordered_map>                         // for unordered_map
#include "boost/graph/detail/adjacency_list.hpp" // for in_degree, out_degree
#include "boost/graph/detail/edge.hpp"           // for operator<, operator!=
#include "boost/graph/named_function_params.hpp" // for bgl_named_params
#include "boost/property_map/property_map.hpp"   // for get, put
#include "boost/tuple/detail/tuple_basic.hpp"    // for get
#include "globalPrimitives.h"                    // for triple, edgeL, stan...
#include "molGraph.h"                            // for atom, molGraph

using namespace std;

/**
 * @brief Halt as soon as any isomorphism is found
 *
 */
struct halting_callback
{
    template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) const
    {
        return false;
    }
};

molGraphBoost edgelistToBoost(molGraph &mg, vector<edgeL> &edgeList, const standardBitset &mask)
{
    molGraphBoost output;
    std::unordered_map<int, int> ht;
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if (mask[i])
        {
            int a = edgeList[i].a, b = edgeList[i].b;
            if (ht.count(a) == 0)
            {
                size_t x = ht.size();
                ht[a] = x;
                string c = mg.mg[a].type;
                add_vertex(atom_vf2(c), output);
            }
            if (ht.count(b) == 0)
            {
                size_t x = ht.size();
                ht[b] = x;
                string c = mg.mg[b].type;
                add_vertex(atom_vf2(c), output);
            }
            int a2 = ht[a], b2 = ht[b];
            add_edge(a2, b2, mg.btype(a, edgeList[i].c), output);
        }
    }
    return output;
}

bool vf2GraphIso(molGraphBoost &mmg, molGraphBoost &tmg)
{
    vertex_comp_t vc =
        boost::make_property_map_equivalent(get(boost::vertex_name, mmg), get(boost::vertex_name, tmg));
    edge_comp_t ec =
        boost::make_property_map_equivalent(get(boost::edge_name, mmg), get(boost::edge_name, tmg));
    halting_callback callback;
    return vf2_graph_iso(mmg, tmg, callback, vertex_order_by_mult(mmg), edges_equivalent(ec).vertices_equivalent(vc));
}