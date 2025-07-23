/**
 * @file molGraph.h
 * @brief code relating to main molecular graph object
 */
#pragma once
#include <stddef.h>           // for size_t
#include <iostream>           // for operator<<, basic_ostream, basic_ostre...
#include <string>             // for allocator, char_traits, basic_string
#include <unordered_map>      // for unordered_map
#include <utility>            // for pair
#include <vector>             // for vector
#include "globalPrimitives.h" // for edgeL, standardBitset, vb

/**
 * @brief Bond struct for molGraph
 */
struct bond
{
    short n;
    short type;
    bond() {}
    bond(short _n, short _type) : n(_n), type(_type) {}
};

/**
 * @brief Atom struct for molGraph
 */
struct atom
{
    std::string type;
    std::vector<bond> list;

    atom() {}
    atom(std::string _type) : type(_type) {}
};

/**
 * @brief Primary graph data structure used in assemblyCpp
 */
struct molGraph
{
    /**
     * @brief Vector of atoms representing nodes in the graph.
     */
    std::vector<atom> mg;
    /**
     * @brief Total number of bonds (edges) in the graph.
     */
    int totalBonds = 0;

    /**
     * @brief Use this function to add atoms/nodes
     * @param _type Type is atom type/node labelling.
     */
    void addAtom(std::string &_type)
    {
        atom a(_type);
        mg.push_back(a);
    }

    /**
     * @brief Use this function to add bonds/edges.
     * @param a Index of first atom/node
     * @param b Index of second atom/node
     * @param type Type is bond order/edge labelling.
     */
    void addBond(int a, int b, short type)
    {
        bond b1(b, type), b2(a, type);
        mg[a].list.push_back(b1);
        mg[b].list.push_back(b2);
        totalBonds++;
    }

    /**
     * @brief Print the graph information to cout
     */
    void printToCout()
    {
        std::cout << "There are " << mg.size() << " atoms in the molecule-graph\n";
        for (size_t i = 0; i < mg.size(); i++)
        {
            std::cout << "Atom " << i + 1 << " is of type " << mg[i].type << " and adjacent to atoms ";
            for (size_t j = 0; j < degree(i); j++)
            {
                std::cout << elem(i, j) + 1 << " with bond order " << btypeS(i, j) << ", ";
            }
            std::cout << '\n';
        }
    }

    /**
     * @brief Get the degree (number of bonds) of atom at index x
     */
    size_t degree(int x)
    {
        return mg[x].list.size();
    }

    short elem(size_t a, size_t b)
    {
        return mg[a].list[b].n;
    }

    /**
     * @brief Get atom type for index i
     */
    std::string atype(size_t i) { return mg[i].type; }

    /**
     * @brief Get bond type as char
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    char btype(size_t a, size_t b)
    {
        if (a < 0 || b < 0)
            return 0;
        return (char)mg[a].list[b].type;
    }

    /**
     * @brief Get bond type as short
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    short btypeS(size_t a, size_t b)
    {
        if (a < 0 || b < 0)
            return 0;
        return mg[a].list[b].type;
    }

    /**
     * @brief Remove all bonds of order 0. Used in preprocessing.
     */
    void collapse()
    {
        size_t originalSize = mg.size(), newSize;
        std::vector<size_t> map, revmap(originalSize, -1);
        molGraph output;
        for (size_t i = 0; i < originalSize; i++)
        {
            revmap[i] = map.size();
            map.push_back(i);
            output.addAtom(mg[i].type);
        }
        newSize = map.size();
        for (size_t i = 0; i < newSize; i++)
        {
            int k = map[i];
            for (size_t j = 0; j < degree(k); j++)
            {
                if (revmap[elem(k, j)] != -1 && btypeS(k, j) != 0)
                {
                    size_t x = revmap[elem(k, j)];
                    if (i < x)
                        output.addBond(i, x, btype(k, j));
                }
            }
        }
        *this = output;
    }

    /**
     * @brief For explicit hydrogen removal
     *
     */
    bool removeAtom(size_t i)
    {
        if (i >= mg.size())
            return false;
        mg[i].type = "COLLAPSE";
        return true;
    }

    /**
     * @brief For explicit hydrogen removal
     *
     */
    void removeAndCollapse()
    {
        size_t originalSize = mg.size(), newSize;
        std::vector<size_t> map, revmap(originalSize, -1);
        molGraph output;
        for (size_t i = 0; i < originalSize; i++)
        {
            if (mg[i].type != "COLLAPSE")
            {
                revmap[i] = map.size();
                map.push_back(i);
                output.addAtom(mg[i].type);
            }
        }
        newSize = map.size();
        for (size_t i = 0; i < newSize; i++)
        {
            int k = map[i];
            for (size_t j = 0; j < degree(k); j++)
            {
                if (revmap[elem(k, j)] != -1 && btypeS(k, j) != 0)
                {
                    size_t x = revmap[elem(k, j)];
                    if (i < x)
                        output.addBond(i, x, btype(k, j));
                }
            }
        }
        *this = output;
    }

    /**
     * @brief Turns molGraph (adjacency list) into equivalent edgelist
     * @return std::vector<edgeL>
     */
    std::vector<edgeL> writeEdgeList()
    {
        std::vector<edgeL> out;
        for (short i = 0; i < mg.size(); i++)
        {
            for (short j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < k)
                {
                    edgeL t(i, k, j);
                    out.push_back(t);
                }
            }
        }
        return out;
    }

    /**
     * @brief For preprocessing, writes edgeList as hash map to detect duplicated bonds
     */
    void writeEdgeList(std::unordered_map<std::string, std::pair<int, edgeL>> &ht)
    {
        for (short i = 0; i < mg.size(); i++)
        {
            for (short j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < k)
                {
                    std::string is = atype(i), ks = atype(k), out;
                    if (is < ks)
                        out = is + btype(i, j) + ks;
                    else
                        out = ks + btype(i, j) + is;
                    edgeL t(i, k, j);
                    if (ht.count(out) == 0)
                    {
                        std::pair<int, edgeL> p(1, t);
                        ht[out] = p;
                    }
                    else
                        ht[out].first++;
                }
            }
        }
    }

    /**
     * @brief Used in preprocessing, removes edges in edgelist
     */
    molGraph negativeEdgeCollapse(std::vector<edgeL> &edgeList)
    {
        for (size_t i = 0; i < edgeList.size(); i++)
        {
            edgeL &el = edgeList[i];
            mg[el.a].list[el.c].type = 0;
        }
        collapse();
        return *this;
    }

    /**
     * @brief For compensating for disjoint fragments in the JAI
     *
     */
    void disjointFragmentsR(vb &visited, int n)
    {
        if (visited[n])
            return;
        visited[n] = 1;
        for (size_t i = 0; i < mg[n].list.size(); i++)
            disjointFragmentsR(visited, elem(n, i));
    }

    /**
     * @brief For compensating for disjoint fragments in the JAI
     *
     * @return int number of disjoint fragments
     */
    int disjointFragments()
    {
        vb visited(mg.size(), 0);
        int count = 0;
        for (size_t i = 0; i < mg.size(); i++)
        {
            if (!visited[i])
            {
                count++;
                disjointFragmentsR(visited, i);
            }
        }
        return count;
    }
};

/**
 * @brief construct new molGraph from input molGraph and boolean edgelist
 *
 * @param mg Input molgraph
 * @param edgeList Corresponding edge list
 * @param mask Input boolean edgelist
 * @param isCyclic Is the graph cyclic, needed for hashing
 * @return molGraph
 */
molGraph constructFromEdgeList(molGraph &mg, std::vector<edgeL> &edgeList,
                               standardBitset &mask, bool &isCyclic);
/**
 * @brief Preprocesses the graph by removing all unique edges for the pathway algorithm
 * @param mg The input molGraph
 * @param writeback The edge list after preprocessing
 * @return molGraph (the final output)
 */
molGraph preprocessWriteback(molGraph &mg, std::vector<edgeL> &writeback);

/// Global variable for the molGraph before and after preprocessing
extern molGraph originalMolecule, targetMolecule;

/**
 * @brief Calls the disjoint-set data structure for the fragmentation function. See Seet et al. section 4.5 for details
 *
 * @param mask Target bitset as input
 * @param maskList List of disjoint bitsets returned
 */
void ufdsMaskConstruct(standardBitset &mask,
                       std::vector<standardBitset> &maskList);