#include "graphHashes.h"
#include <algorithm>          // for min
#include <bitset>             // for hash
#include <string>             // for basic_string, operator==, hash, string
#include <unordered_map>      // for unordered_map
#include <utility>            // for pair
#include <vector>             // for vector
#include "globalPrimitives.h" // for atypeHash, bitsetHashTable, standardBi...
#include "molGraph.h"         // for molGraph, atom, constructFromEdgeList
#include "treeCanon.h"        // for centroidTreeCanon

using namespace std;

graphHash::graphHash(molGraph &mg, int depth, bool isCyclic, standardBitset &_mask)
{
    mask = _mask;
    if (isCyclic)
    {
        hashes.resize(mg.mg.size(), 0);
        calcHash(mg, depth);
    }
    else
        treeHash = centroidTreeCanon(mg, 0);
}

void graphHash::calcHash(molGraph &mg, int _depth)
{
    const double depthFactor = 0.33;
    int depth = min(_depth, HASH_DEPTH_MAX);
    vector<double> dhashes(hashes.size(), 0.0);
    for (size_t i = 0; i < mg.mg.size(); i++)
    {
        string atype1 = mg.mg[i].type;
        if (atypeHash.count(atype1) == 0)
            atypeHash[atype1] = (atypeHash.size() + 1) * 5;
        dhashes[i] = (0.0 + atypeHash[atype1]) / depthFactor;
        for (size_t j = 0; j < mg.degree(i); j++)
        {
            short btype = mg.btype(i, j);
            string atype = mg.atype(mg.elem(i, j));
            if (atypeHash.count(atype) == 0)
                atypeHash[atype] = (atypeHash.size() + 1) * 5;
            int atypeInt = atypeHash[atype];
            dhashes[i] += atypeInt + btype;
        }
    }
    for (int k = 1; k < depth; k++)
    {
        vector<double> oldHashes = dhashes;
        for (size_t i = 0; i < mg.mg.size(); i++)
        {
            for (size_t j = 0; j < mg.degree(i); j++)
            {
                short btype = mg.btype(i, j);
                dhashes[i] += (oldHashes[mg.elem(i, j)] + btype) * k * depthFactor;
            }
        }
    }
    for (size_t i = 0; i < hashes.size(); i++)
        hashes[i] = (float)dhashes[i];
}

std::unordered_map<graphHash, pii> graphHashMap;

int canonise(standardBitset &mask)
{
    bool isCyclic;
    vector<edgeL> &edgeList = univEdgeList;
    size_t x;
    if (bitsetHashTable.count(mask) == 0)
    {
        molGraph mg = constructFromEdgeList(targetMolecule, edgeList, mask, isCyclic);
        graphHash g(mg, mg.mg.size(), isCyclic, mask);
        if (graphHashMap.count(g) == 0)
        {
            x = graphHashMap.size();
            graphHashMap[g].first = x;
            graphHashMap[g].second = 1;
        }
        else
        {
            x = graphHashMap[g].first;
            graphHashMap[g].second++;
        }
        bitsetHashTable[mask].first = x;
        bitsetHashTable[mask].second = graphHashMap[g].second;
    }
    else
    {
        x = bitsetHashTable[mask].first;
    }
    return x;
}