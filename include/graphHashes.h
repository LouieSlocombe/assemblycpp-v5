
/**
 * @file graphHashes.h
 * @brief code relating to hashing graph objects
 */

#pragma once
#include <algorithm>          // for sort
#include <cstddef>            // for size_t
#include <string>             // for hash, operator==, string, __str_hash_base
#include <unordered_map>      // for hash, unordered_map
#include <variant>            // for hash
#include <vector>             // for vector
#include "globalPrimitives.h" // for standardBitset, univEdgeList, pii
#include "molGraph.h"         // for molGraph (ptr only), targetMolecule
#include "vf2.h"              // for edgelistToBoost, vf2GraphIso, molGraph...

/**
 * @brief Hashes a molecular graph
 */
struct graphHash
{
    /// graph to hash expressed as a bitset of edges
    standardBitset mask;
    /// hashes stored as floats
    std::vector<float> hashes;
    /// if the graph is acyclic, the tree hash function is used, and the output stored here
    std::string treeHash;

    graphHash() {}

    /**
     * @brief Construct a new graph Hash object
     *
     * @param mg molGraph to be hashed
     * @param depth Depth of the BFS
     * @param isCyclic Is the molecule cyclic
     * @param _mask Boolean edgelist of the molGraph
     */
    graphHash(molGraph &mg, int depth, bool isCyclic, standardBitset &_mask);

    /**
     * @brief Calculate hash for subgraph unordered_map using BFS approach
     *
     * @param mg molGraph to be hashed
     * @param _depth Depth of the BFS
     */
    void calcHash(molGraph &mg, int _depth);

    /**
     * @brief Check isomorphism between two graphs. Uses tree isomorphism if acyclic, else uses vf2 subgraph isomorphism
     *
     * @param g2 other graph to be compared
     * @return true if graphs are isomorphic
     * @return false otherwise
     */
    bool operator==(const graphHash &g2) const
    {
        if (treeHash.length() != g2.treeHash.length())
            return false;
        if (treeHash.length() == 0)
        {
            molGraphBoost g1mg = edgelistToBoost(targetMolecule, univEdgeList, this->mask),
                          g2mg = edgelistToBoost(targetMolecule, univEdgeList, g2.mask);
            return vf2GraphIso(g1mg, g2mg);
        }
        else
        {
            return treeHash == g2.treeHash;
        }
    }
};

/**
 * @brief Hash for subgraph unordered_map
 */
template <>
struct std::hash<graphHash>
{
    size_t operator()(const graphHash &gh) const
    {
        if (gh.treeHash.length() == 0)
        {
            std::vector<float> sortedHashes = gh.hashes;
            sort(sortedHashes.begin(), sortedHashes.end());
            size_t res = 17;
            for (size_t i = 0; i < sortedHashes.size(); i++)
            {
                int k = int(sortedHashes[i] * 1024);
                res = res * 31 + std::hash<int>()(k);
            }
            return res;
        }
        else
        {
            std::hash<string> hasher;
            return hasher(gh.treeHash);
        }
    }
};

/**
 * @brief TODO: document
 */
extern std::unordered_map<graphHash, pii> graphHashMap;

/**
 * @brief Returns unique hash val for subgraph. See Seet et al. section 4.3 Enumeration
 *
 * @param mask Boolean edgelist to be canonised
 * @return int canonical value
 */
int canonise(standardBitset &mask);