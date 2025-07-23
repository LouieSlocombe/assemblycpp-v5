#include "fragmentation.h"
#include <bitset>              // for bitset
#include <cstddef>             // for size_t, std
#include <unordered_map>       // for unordered_map
#include <unordered_set>       // for unordered_set, operator!=, _Node_iter...
#include <vector>              // for vector
#include "assemblyState.h"     // for assemblyState, apWrapper, pathAssembl...
#include "duplicateMatching.h" // for validMatchings
#include "globalPrimitives.h"  // for standardBitset, bitsetHashTable
#include "graphHashes.h"       // for canonise
#include "molGraph.h"          // for ufdsMaskConstruct

using namespace std;

void fragmentAssemblyState(assemblyState &_target, validMatchings &matching,
                           assemblyState &_result)
{
    vector<standardBitset> &masks = _target.masks;
    standardBitset f1 = matching.first, f2 = matching.second;
    bool same = 1;
    if (matching.frag1 != matching.frag2)
        same = 0;
    _result.masks.push_back(f1);
    if (same)
    {
        standardBitset resultMask = masks[matching.frag1];
        resultMask ^= f1;
        resultMask ^= f2;
        ufdsMaskConstruct(resultMask, _result.masks);
    }
    else
    {
        standardBitset resultMask1 = masks[matching.frag1];
        resultMask1 ^= f1;
        ufdsMaskConstruct(resultMask1, _result.masks);
        standardBitset resultMask2 = masks[matching.frag2];
        resultMask2 ^= f2;
        ufdsMaskConstruct(resultMask2, _result.masks);
    }
    for (size_t i = 0; i < _result.masks.size(); i++)
    {
        canonise(_result.masks[i]);
    }
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (i != matching.frag1 && i != matching.frag2 && masks[i] != 0)
        {
            vector<standardBitset> tempMasks;
            if (bitsetHashTable.count(masks[i]) == 0)
            {
                ufdsMaskConstruct(masks[i], tempMasks);
                for (size_t j = 0; j < tempMasks.size(); j++)
                {
                    canonise(tempMasks[j]);
                    _result.masks.push_back(tempMasks[j]);
                }
            }
            else
                _result.masks.push_back(masks[i]);
        }
    }
}

void clearPathMap()
{
    vector<std::unordered_set<assemblyPath *>::iterator> vit;
    for (auto it = pathAssemblyMap.begin(); it != pathAssemblyMap.end(); ++it)
    {
        delete (*it).ap;
    }
    pathAssemblyMap.clear();
}