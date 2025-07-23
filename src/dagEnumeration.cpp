#include "dagEnumeration.h"
#include <algorithm>          // for sort
#include <cstddef>            // for size_t, std
#include <unordered_map>      // for unordered_map, _Node_iterator, operator!=
#include <utility>            // for pair
#include <vector>             // for vector
#include "globalPrimitives.h" // for standardBitset, bitsetHashTable

using namespace std;

vector<vector<dagNode>> DAG;

void convertDag(vector<std::unordered_map<standardBitset, pair<int, vector<standardBitset>>>> &tempDag)
{
    DAG.resize(tempDag.size());
    vector<std::unordered_map<standardBitset, int>> bitsetToIndex(tempDag.size());
    for (size_t i = 0; i < tempDag.size() - 1; i++)
    {
        for (auto it = tempDag[i].begin(); it != tempDag[i].end(); ++it)
        {
            vector<standardBitset> &list = it->second.second;
            size_t trueSize = list.size();
            for (size_t j = 0; j < list.size(); j++)
            {
                std::unordered_map<standardBitset, pair<int, vector<standardBitset>>> &nextMap = tempDag[i + 1];
                if (nextMap.count(list[j]) == 0)
                {
                    list[j] = 0;
                    trueSize--;
                }
            }
            vector<standardBitset> trueList(trueSize);
            size_t k = 0;
            for (size_t j = 0; j < list.size(); j++)
            {
                if (list[j] != 0)
                {
                    trueList[k] = list[j];
                    k++;
                }
            }
            it->second.second = trueList;
            if (i > 0)
            {
                it->second.first = bitsetHashTable[it->first].first;
                size_t x = bitsetToIndex[i].size();
                bitsetToIndex[i][it->first] = x;
            }
        }
    }
    for (auto it = tempDag[tempDag.size() - 1].begin(); it != tempDag[tempDag.size() - 1].end(); ++it)
    {
        it->second.first = bitsetHashTable[it->first].first;
    }
    for (size_t i = 0; i < DAG.size() - 1; i++)
    {
        for (auto it = tempDag[i].begin(); it != tempDag[i].end(); ++it)
        {
            dagNode dn(it->first, it->second.first, it->second.second, bitsetToIndex[i + 1]);
            DAG[i].push_back(dn);
        }
    }
    sort(DAG[0].begin(), DAG[0].end(), CompareDagNode());
}
