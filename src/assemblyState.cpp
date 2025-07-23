#include "assemblyState.h"
#include <stddef.h>      // for size_t
#include <algorithm>     // for max, sort
#include <bitset>        // for bitset, operator<<, hash
#include <cmath>         // for ceil, log2
#include <fstream>       // for operator<<, basic_ostream, basic_ostream::o...
#include <iostream>      // for cout
#include <locale>        // for num_get, num_put, numpunct
#include <string>        // for char_traits
#include <unordered_map> // for unordered_map
#include <unordered_set> // for unordered_set
#include <utility>       // for pair
#include <vector>        // for vector
#include "molGraph.h"    // for constructFromEdgeList, molGraph, targetMole...

using namespace std;

assemblyPath *minAssemblyPath = nullptr;

std::unordered_set<apWrapper> pathAssemblyMap;

int assemblyState::maxFragSizeF()
{
    return masks[0].count();
}

int assemblyState::maxDupBonds()
{
    int dupBonds2 = 0, dupBondsTotal, maxFragSize = maxFragSizeF();
    vi sizeList(masks.size());

    for (size_t i = 0; i < masks.size(); i++)
    {
        sizeList[i] = masks[i].count();
    }

    for (size_t i = 0; i < sizeList.size(); i++)
        dupBonds2 += sizeList[i] / 2;
    dupBonds2--;
    for (int j = 3; j <= maxFragSize; j++)
    {
        dupBondsTotal = 0;
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            dupBondsTotal += (sizeList[i] - sizeList[i] / j);
            if (sizeList[i] % j != 0)
                dupBondsTotal--;
        }
        dupBondsTotal -= ceil(log2(j));
        if (dupBondsTotal > dupBonds2)
            dupBonds2 = dupBondsTotal;
    }
    return dupBonds2;
}

int assemblyState::maxDupBonds(vi &sizeListMain, int maxFragSize, vi &sizeList)
{
    int dupBonds2 = 0, dupBondsTotal;

    int j = maxFragSize;
    vi adjustedSizeList(sizeList.size()), adjustedSizeList2(sizeList.size());
    for (size_t i = 0; i < sizeList.size(); i++)
    {
        adjustedSizeList[i] = sizeList[i] - sizeList[i] % j;
        adjustedSizeList2[i] = sizeListMain[i] - adjustedSizeList[i];
    }
    dupBondsTotal = 0;
    for (size_t i = 0; i < sizeList.size(); i++)
    {
        dupBondsTotal += (adjustedSizeList[i] - adjustedSizeList[i] / j);
        dupBondsTotal += (adjustedSizeList2[i] - adjustedSizeList2[i] / (j - 1));
        if (adjustedSizeList2[i] % (j - 1) != 0)
            dupBondsTotal--;
    }
    dupBondsTotal -= ceil(log2(j));
    return dupBondsTotal;
}

int assemblyState::maxDupBonds(vi &sizeListMain, int maxFragSize, vector<standardBitset> &targetMasks)
{
    int dupBonds2 = 0, dupBondsTotal;
    vi sizeList(targetMasks.size());

    for (size_t i = 0; i < masks.size(); i++)
    {
        sizeList[i] = targetMasks[i].count();
    }

    int j = maxFragSize;
    vi adjustedSizeList(sizeList.size()), adjustedSizeList2(sizeList.size());
    for (size_t i = 0; i < sizeList.size(); i++)
    {
        adjustedSizeList[i] = sizeList[i] - sizeList[i] % j;
        adjustedSizeList2[i] = sizeListMain[i] - adjustedSizeList[i];
    }
    dupBondsTotal = 0;
    for (size_t i = 0; i < sizeList.size(); i++)
    {
        dupBondsTotal += (adjustedSizeList[i] - adjustedSizeList[i] / j);
        dupBondsTotal += (adjustedSizeList2[i] - adjustedSizeList2[i] / (j - 1));
        if (adjustedSizeList2[i] % (j - 1) != 0)
            dupBondsTotal--;
    }
    dupBondsTotal -= ceil(log2(j));
    return dupBondsTotal;
}

void assemblyState::maxDupBonds(vi &fragSizeList, int maxFragSize, vector<vector<standardBitset>> &targetMasks)
{
    int dupBonds2 = 0, dupBondsTotal;
    fragSizeList.resize(maxFragSize - 1);
    vector<vi> sizeLists(fragSizeList.size());

    for (size_t j = 0; j < targetMasks.size(); j++)
    {
        sizeLists[j].assign(targetMasks[j].size(), 0);
        for (size_t i = 0; i < masks.size(); i++)
        {
            sizeLists[j][i] = targetMasks[j][i].count();
        }
    }

    for (size_t i = 0; i < sizeLists[0].size(); i++)
        dupBonds2 += sizeLists[0][i] / 2;
    dupBonds2--;
    fragSizeList[0] = dupBonds2;

    for (int j = 3; j <= maxFragSize; j++)
    {
        vi &sizeList = sizeLists[j - 2], &sizeList2 = sizeLists[0];
        vi adjustedSizeList(sizeList.size()), adjustedSizeList2(sizeList.size());
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            adjustedSizeList[i] = sizeList[i] - sizeList[i] % j;
            adjustedSizeList2[i] = sizeList2[i] - adjustedSizeList[i];
        }
        dupBondsTotal = 0;
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            dupBondsTotal += (adjustedSizeList[i] - adjustedSizeList[i] / j);
            dupBondsTotal += (adjustedSizeList2[i] - adjustedSizeList2[i] / (j - 1));
            if (adjustedSizeList2[i] % (j - 1) != 0)
                dupBondsTotal--;
        }
        dupBondsTotal -= ceil(log2(j));
        fragSizeList[j - 2] = dupBondsTotal;
    }
}

int assemblyState::maxDupBonds(int maxFragSize)
{
    int dupBonds2 = 0, dupBondsTotal;
    vi sizeList(masks.size());

    for (size_t i = 0; i < masks.size(); i++)
    {
        sizeList[i] = masks[i].count();
    }

    for (size_t i = 0; i < sizeList.size(); i++)
        dupBonds2 += sizeList[i] / 2;
    dupBonds2--;
    for (int j = 3; j <= maxFragSize; j++)
    {
        dupBondsTotal = 0;
        for (size_t i = 0; i < sizeList.size(); i++)
        {
            dupBondsTotal += (sizeList[i] - sizeList[i] / j);
            if (sizeList[i] % j != 0)
                dupBondsTotal--;
        }
        dupBondsTotal -= ceil(log2(j));
        if (dupBondsTotal > dupBonds2)
            dupBonds2 = dupBondsTotal;
    }
    return dupBonds2;
}

int assemblyState::lowBoundAI()
{
    return totalBonds - sumDupBonds - 1 - maxDupBonds();
}

int assemblyState::lowBoundAI(int maxFragSize, int estimate)
{
    if (maxFragSize > 2)
        estimate = max(estimate, maxDupBonds(maxFragSize - 1));
    return totalBonds - sumDupBonds - 1 - estimate;
}

int assemblyState::AI()
{
    return totalBonds - sumDupBonds - 1;
}

vi assemblyState::assemblyHashCalculator()
{
    vi sorted(masks.size(), -1);
    for (size_t i = 0; i < masks.size(); i++)
    {
        if (bitsetHashTable.count(masks[i]))
            sorted[i] = bitsetHashTable[masks[i]].first;
    }
    sort(sorted.begin() + 1, sorted.end());
    return sorted;
}

void assemblyState::print()
{
    cout << "printing assembly state: " << ix << '\n';
    cout << "masks:\n";
    for (size_t i = 0; i < masks.size(); i++)
    {
        cout << masks[i] << '\n';
    }
    for (size_t i = 0; i < masks.size(); i++)
    {
        cout << "Fragment: " << i << '\n';
        bool isCyclic;
        molGraph mg = constructFromEdgeList(targetMolecule, univEdgeList, masks[i], isCyclic);
        mg.printToCout();
    }
    cout << "AI: " << AI() << '\n';
    cout << "lowbound: " << lowBoundAI() << '\n';
}
