#include <string>        // std::string
#include <vector>        // std::vector
#include <iostream>      // std::cout
#include <unordered_map> // std::unordered_map
#include "molGraph.h"
#include "globalPrimitives.h" // standardBitset, edgeL, originalMolecule, targetMolecule, univEdgeList
#include "ufds.h"             // disjointSet, ufdsSplit

using namespace std;

molGraph constructFromEdgeList(molGraph &mg, vector<edgeL> &edgeList,
                               standardBitset &mask, bool &isCyclic)
{
    disjointSet u(mg.mg.size());
    molGraph output;
    std::unordered_map<int, int> ht;
    isCyclic = 0;
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if (mask[i] != 0)
        {
            int a = edgeList[i].a, b = edgeList[i].b, c = 0;
            if (ht.count(a) == 0 && ht.count(b) == 0)
            {
                size_t x = ht.size();
                ht[a] = x;
                output.addAtom(mg.mg[a].type);
                size_t y = ht.size();
                ht[b] = y;
                output.addAtom(mg.mg[b].type);
                u.insert(x, x);
                u.insert(y, x);
            }
            else
            {
                if (ht.count(a) == 0)
                {
                    size_t x = ht.size();
                    ht[a] = x;
                    output.addAtom(mg.mg[a].type);
                    u.insert(x, ht[b]);
                }
                else
                    c++;
                if (ht.count(b) == 0)
                {
                    size_t x = ht.size();
                    ht[b] = x;
                    output.addAtom(mg.mg[b].type);
                    u.insert(x, ht[a]);
                }
                else
                    c++;
            }
            int a2 = ht[a], b2 = ht[b];
            output.addBond(a2, b2, mg.btype(a, edgeList[i].c));
            if (c == 2)
            {
                isCyclic |= u.merge(a2, b2);
            }
        }
    }
    return output;
}

molGraph preprocessWriteback(molGraph &mg, vector<edgeL> &writeback)
{
    std::unordered_map<string, pair<int, edgeL>> ht;
    molGraph out = mg;
    mg.writeEdgeList(ht);
    vector<edgeL> v;
    for (auto it = ht.begin(); it != ht.end(); ++it)
    {
        if (it->second.first == 1)
        {
            v.push_back(it->second.second);
        }
    }
    out.negativeEdgeCollapse(v);
    writeback = v;
    return out;
}

/// Global variable for the molGraph before and after preprocessing
molGraph originalMolecule, targetMolecule;

void ufdsMaskConstruct(standardBitset &mask,
                       vector<standardBitset> &maskList)
{
    vector<edgeL> &edgeList = univEdgeList;
    ufdsSplit u(targetMolecule.mg.size());
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if (mask[i] != 0)
        {
            int a = edgeList[i].a, b = edgeList[i].b;
            if (u.elements[a].parent == -1 && u.elements[b].parent == -1)
            {
                u.doubleInsert(b, a, i);
            }
            else
            {
                if (u.elements[a].parent == -1)
                {
                    u.insert(a, b, i);
                }
                else if (u.elements[b].parent == -1)
                {
                    u.insert(b, a, i);
                }
                else
                {
                    u.merge(a, b, i);
                }
            }
        }
    }
    u.split(maskList);
}