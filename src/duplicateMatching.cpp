#include "duplicateMatching.h"
#include <bitset>             // for bitset, operator&, operator|, hash
#include <cstddef>            // for size_t, std
#include <map>                // for map, operator==, _Rb_tree_iterator
#include <unordered_map>      // for unordered_map
#include <unordered_set>      // for unordered_set
#include <utility>            // for pair
#include <vector>             // for vector
#include "dagEnumeration.h"   // for dagNode, DAG
#include "globalPrimitives.h" // for standardBitset, triple, univEdgeList

using namespace std;

initialPotentialDuplicate::initialPotentialDuplicate(int x, standardBitset &_fragMask, size_t _fragment)
{
    fragMask = _fragMask;
    fragment = _fragment;
    mask.set(x);
    atomMask.set(univEdgeList[x].a);
    atomMask.set(univEdgeList[x].b);
}

void initialPotentialDuplicate::generate(vector<initialPotentialDuplicate> &q, size_t fragment, std::unordered_set<standardBitset> &maskMap)
{
    vector<edgeL> &edgeList = univEdgeList;
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if ((mask[i] == 0) && (fragMask[i] != 0))
        {
            standardBitset temp1 = 0, temp2 = 0;
            temp1.set(edgeList[i].a);
            temp2.set(edgeList[i].b);
            if ((temp1 & atomMask) == temp1 || (temp2 & atomMask) == temp2)
            {
                standardBitset tempMask = mask;
                tempMask.set(i);
                if (maskMap.count(tempMask) == 0)
                {
                    maskMap.insert(tempMask);
                    initialPotentialDuplicate g = *this;
                    g.mask.set(i);
                    g.atomMask |= (temp1 | temp2);
                    q.push_back(g);
                }
            }
        }
    }
}

void initialPotentialDuplicate::generateDAG(vector<initialPotentialDuplicate> &q, size_t fragment, std::unordered_set<standardBitset> &maskMap,
                                            vector<std::unordered_map<standardBitset, pair<int, vector<standardBitset>>>> &tempDag)
{
    vector<edgeL> &edgeList = univEdgeList;
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if ((mask[i] == 0) && (fragMask[i] != 0))
        {
            standardBitset temp1, temp2;
            temp1.set(edgeList[i].a),
                temp2.set(edgeList[i].b);
            if ((temp1 & atomMask) == temp1 || (temp2 & atomMask) == temp2)
            {
                standardBitset tempMask = mask;
                tempMask.set(i);
                if (maskMap.count(tempMask) == 0)
                {
                    maskMap.insert(tempMask);
                    initialPotentialDuplicate g = *this;
                    g.mask.set(i);
                    g.atomMask |= (temp1 | temp2);
                    q.push_back(g);
                    vector<standardBitset> &adjList = tempDag[mask.count() - 1][mask].second;
                    pair<int, vector<standardBitset>> p;
                    tempDag[mask.count()][g.mask] = p;
                    adjList.push_back(g.mask);
                }
            }
        }
    }
}

bool dagGenerate(potentialDuplicate &d, map<int, dagDuplicateSet> &stmap, standardBitset &fragment,
                 size_t size, int ordinal, size_t frags)
{
    bool overweight = 0;
    for (size_t i = 0; i < DAG[size - 1][d.idx].children.size(); i++)
    {
        dagNode &dn = DAG[size][DAG[size - 1][d.idx].children[i]];
        if (((dn.mask | fragment) == fragment))
        {
            if (dn.ix <= ordinal)
            {
                auto it = stmap.find(dn.ix);
                if (it == stmap.end())
                {
                    dagDuplicateSet ss(size + 1, frags);
                    ss.insert(potentialDuplicate(dn.mask, d.fragment, DAG[size - 1][d.idx].children[i]));
                    stmap[dn.ix] = ss;
                }
                else
                {
                    it->second.insert(potentialDuplicate(dn.mask, d.fragment, DAG[size - 1][d.idx].children[i]));
                }
            }
            else
                overweight = 1;
        }
    }
    return overweight;
}

bool dagDuplicateGenerator(dagDuplicateSet &ds, map<int, dagDuplicateSet> &stmap,
                           vector<standardBitset> &takenMasks, vector<standardBitset> &stateMasks, int ordinal, bool &overweight, bool last)
{
    bool output = 0;
    vb alive(ds.list.size(), 0);
    for (size_t i = 0; i < ds.list.size(); i++)
    {
        size_t frag = ds.list[i].fragment;
        if (ds.size > 0)
        {
            for (size_t j = i + 1; j < ds.list.size(); j++)
            {
                if (frag == ds.list[j].fragment)
                {
                    if ((ds.list[i].mask & ds.list[j].mask) == 0)
                    {
                        alive[i] = 1;
                        alive[j] = 1;
                    }
                }
                else
                {
                    alive[i] = 1;
                    alive[j] = 1;
                }
            }
        }
        else
            alive[i] = 1;
        if (alive[i])
        {
            takenMasks[frag] |= ds.list[i].mask;
            ds.dead = 0;
            if (!last)
            {
                overweight |= dagGenerate(ds.list[i], stmap, stateMasks[frag], ds.size,
                                          ordinal, ds.maskList.size());
                output = 1;
            }
        }
    }
    return output;
}